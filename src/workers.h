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

#ifndef _workers_h
#define _workers_h
	
	typedef pthread_mutex_t workers_mutex;  /* Redefinition of the mutex. */
	
	typedef struct{
		pthread_t *threads;           /* Array of workers created.   */
		pthread_mutex_t excluder;     /* Allows to control the pool. */
		pthread_cond_t newtask;       /* Sleeps the idle workers.    */
		pthread_cond_t launcher;      /* Control launcher function.  */
		pthread_cond_t waiter;        /* Control barrier function.   */
		number num,stop,next,exec;    /* Control variables of works. */
		void (*routine)(void*),*data; /* The next work to execute.   */
	}workers;
	
	/* Public function prototypes to work with the workers pool. */
	workers *workers_create(number);
	void workers_addtask(workers*,void(*)(void*),void*);
	void workers_waitall(workers*);
	void workers_finish(workers*);
	number workers_number(workers*);
	
	/* Public function prototypes to work with the workers mutex. */
	workers_mutex *workers_mutex_create();
	void workers_mutex_delete(workers_mutex*);
	void workers_mutex_lock(workers_mutex*);
	void workers_mutex_unlock(workers_mutex*);
	
#endif
