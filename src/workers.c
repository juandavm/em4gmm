/* Expectation Maximization for Gaussian Mixture Models.
Copyright (C) 2012-2022 Juan Valor

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#include "global.h"
#include "workers.h"

/* This are the worker threads that executes the works. */
void *workers_thread(void *tdata){
	workers *t=(workers*)tdata;
	void (*routine_exec)(void*),*data;
	pthread_mutex_lock(&t->excluder);
	while(1){
		if(t->next==1){
			routine_exec=t->routine,data=t->data,t->next=0;
			pthread_mutex_unlock(&t->excluder);
			pthread_cond_signal(&t->launcher);
			routine_exec(data);
			pthread_mutex_lock(&t->excluder);
			if((--t->exec)==0)
				pthread_cond_signal(&t->waiter);
		}else if(t->stop!=1)
			pthread_cond_wait(&t->newtask,&t->excluder);
		else if(t->stop==1)break;
	}
	pthread_mutex_unlock(&t->excluder);
	return NULL;
}

/* Return the number of workers on the pool. */
number workers_number(workers *pool){
	return pool->num;
}

/* A mutex lock for the worker that uses it. */
void workers_mutex_lock(workers_mutex *mutex){
	pthread_mutex_lock(mutex);
}

/* A mutex unlock for the worker that uses it. */
void workers_mutex_unlock(workers_mutex *mutex){
	pthread_mutex_unlock(mutex);
}

/* Creates and initializes a mutex, and returns it. */
workers_mutex *workers_mutex_create(){
	workers_mutex *mutex=(workers_mutex*)calloc(1,sizeof(workers_mutex));
	pthread_mutex_init(mutex,NULL);
	return mutex;
}

/* Destroy and free the mutex passed as argument. */
void workers_mutex_delete(workers_mutex *mutex){
	pthread_mutex_destroy(mutex);
	free(mutex);
}

/* Barrier that waits the end of all the workers. */
void workers_waitall(workers *pool){
	pthread_mutex_lock(&pool->excluder);
	if(pool->exec!=0)
		pthread_cond_wait(&pool->waiter,&pool->excluder);
	pthread_mutex_unlock(&pool->excluder);
}

/* Adds one task to the workers, or waits if it cannot. */
void workers_addtask(workers *pool,void (*routine)(void*),void *data){
	pthread_mutex_lock(&pool->excluder);
	while(1){
		if(pool->next==0){
			pool->data=data,pool->routine=routine;
			pool->next=1,pool->exec++;
			pthread_mutex_unlock(&pool->excluder);
			pthread_cond_signal(&pool->newtask);
			break;
		}else pthread_cond_wait(&pool->launcher,&pool->excluder);
	}
}

/* Creates a pool of worker threads to compute faster. */
workers *workers_create(number num){
	workers *pool=(workers*)calloc(1,sizeof(workers));
	pool->threads=(pthread_t*)calloc(pool->num=num,sizeof(pthread_t));
	pthread_mutex_init(&pool->excluder,NULL);
	pthread_cond_init(&pool->newtask,NULL);
	pthread_cond_init(&pool->launcher,NULL);
	pthread_cond_init(&pool->waiter,NULL);
	number i; for(i=0;i<num;i++)
		pthread_create(&pool->threads[i],NULL,workers_thread,(void*)pool);
	return pool;
}

/* Ends the execution and free the resources. */
void workers_finish(workers *pool){
	pthread_mutex_lock(&pool->excluder);
	number i; pool->stop=1;
	pthread_mutex_unlock(&pool->excluder);
	pthread_cond_broadcast(&pool->newtask);
	for(i=0;i<pool->num;i++)
		pthread_join(pool->threads[i],NULL);
	pthread_cond_destroy(&pool->newtask);
	pthread_cond_destroy(&pool->launcher);
	pthread_cond_destroy(&pool->waiter);
	pthread_mutex_destroy(&pool->excluder);
	free(pool->threads);
	free(pool);
}
