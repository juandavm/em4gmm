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
#include "data.h"

#define BUFFER_SIZE 2*1024*1024 /* Loader samples cache size.   */

typedef struct{
	data *feas;            /* Shared pointer to samples structure.  */
	number r,s,d,c,sign;   /* Internal variables shared by threads. */
	number header,point;   /* Internal variables shared by threads. */
	decimal *aux,next,dec; /* Internal variables shared by threads. */
	char *buff;            /* Pointer to loaded bytes on memory.    */
}loader;

/* Asynchronous implementation of the buffer to feas converter. */
void thread_loader(void *tdata){
	loader *t=(loader*)tdata; number i;
	for(i=0;i<t->r;i++){
		if(t->buff[i]>='0'){ /* Convert an array into a decimal number. */
			t->next=(t->point==0)?((10*t->next)+t->buff[i]-'0'):
				(t->next+((t->dec*=0.1)*(t->buff[i]-'0'))),t->c++;
		}else if(t->buff[i]<'-'&&t->c>0){
			if(t->header==0){
				if(t->d==0){ /* If we start a new sample, set the memory. */
					t->feas->data[t->s]=t->aux;
					t->aux+=t->feas->dimension;
				}
				t->feas->data[t->s][t->d]=t->next*t->sign;
				t->feas->mean[t->d]+=t->feas->data[t->s][t->d];
				t->feas->variance[t->d]+=t->feas->data[t->s][t->d]*t->feas->data[t->s][t->d];
				if((++t->d)==t->feas->dimension)t->d=0,t->s++;
			}else if(t->header==1){ /* Finish header reading and alloc needed memory. */
				t->feas->samples=(int)t->next;
				t->feas->mean=(decimal*)calloc(2*t->feas->dimension,sizeof(decimal));
				t->feas->variance=t->feas->mean+t->feas->dimension;
				t->feas->data=(decimal**)calloc(t->feas->samples,sizeof(decimal*));
				t->aux=(decimal*)calloc(t->feas->dimension*t->feas->samples,sizeof(decimal));
				t->header=t->s=t->d=0;
			}else if(t->header==2)t->feas->dimension=(int)t->next,t->header=1;
			t->sign=t->dec=1,t->next=t->point=t->c=0;
		}else if(t->buff[i]=='-')t->sign=-1;
		else if(t->buff[i]=='.')t->point=1;
	}
}

/* Load the samples from a plain text file with the specified format. */
data *feas_load(char *filename,workers *pool){
	data *feas=(data*)calloc(1,sizeof(data));
	loader *t=(loader*)calloc(1,sizeof(loader));
	t->s=t->d=t->c=t->point=t->next=0,t->header=2,t->dec=t->sign=1;
	t->feas=feas,t->buff=NULL; number i,r;
	gzFile f=gzopen(filename,"rb"); /* Read the file using zlib library. */
	if(!f)fprintf(stderr,"Error: Not %s feature file found.\n",filename),exit(1);
	while(!gzeof(f)){
		char *buff=(char*)calloc(BUFFER_SIZE,sizeof(char));
		r=gzread(f,buff,BUFFER_SIZE); /* Read the buffer and do asynchronous load. */
		if(t->buff!=NULL){
			workers_waitall(pool);
			free(t->buff);
		}
		t->buff=buff,t->r=r;
		workers_addtask(pool,thread_loader,(void*)t);
	}
	workers_waitall(pool);
	gzclose(f); free(t->buff); free(t);
	for(i=0;i<feas->dimension;i++){ /* Compute the mean and variance of the data. */
		feas->mean[i]/=feas->samples;
		feas->variance[i]=(feas->variance[i]/feas->samples)-(feas->mean[i]*feas->mean[i]);
		if(feas->variance[i]<0.000001)feas->variance[i]=0.000001;
	}
	return feas;
}

/* Free the allocated memory of the samples. */
void feas_delete(data *feas){
	free(feas->data[0]);
	free(feas->data);
	free(feas->mean);
	free(feas);
}
