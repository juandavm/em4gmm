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
#include "gmm.h"

typedef struct{
	workers_mutex *mutex; /* Common mutex to lock shared data.   */
	data *feas;           /* Shared pointer to loaded samples.   */
	gmm *gmix;            /* Shared pointer to gaussian mixture. */
	number ini, end;      /* Initial and final sample processed. */
}trainer;

/* Parallel implementation of the E Step of the EM algorithm. */
void thread_trainer(void *tdata){
	trainer *t=(trainer*)tdata; /* Get the data for the thread and alloc memory. */
	gmm *gmix=t->gmix; data *feas=t->feas;
	decimal *zval=(decimal*)calloc(2*gmix->num,sizeof(decimal)),*prob=zval+gmix->num;
	decimal *mean=(decimal*)calloc(2*gmix->num*gmix->dimension,sizeof(decimal)),llh=0;
	decimal *dcov=mean+(gmix->num*gmix->dimension),x,tz,rmean,max,mep=log(DBL_EPSILON);
	number *tfreq=(number*)calloc(gmix->num,sizeof(number)),i,j,m,inc,c;
	for(i=t->ini;i<t->end;i++){
		max=-HUGE_VAL,c=-1;
		for(m=0;m<gmix->num;m++){ /* Compute expected class value of the sample. */
			prob[m]=gmix->mix[m].cgauss;
			for(j=0;j<gmix->dimension;j++){
				x=feas->data[i][j]-gmix->mix[m].mean[j];
				prob[m]-=(x*x)*gmix->mix[m].dcov[j];
			}
			prob[m]*=0.5;
			if(max<prob[m])max=prob[m],c=m;
		}
		for(m=0,x=0;m<gmix->num;m++){ /* Do not use Viterbi aproximation. */
			rmean=prob[m]-max;
			if(rmean>mep)x+=exp(rmean); /* Use machine epsilon to avoid make exp's. */
		}
		llh+=(rmean=max+log(x)); tfreq[c]++;
		for(m=0;m<gmix->num;m++){ /* Accumulate counts of the sample in memory. */
			if((x=prob[m]-rmean)>mep){ /* Use machine epsilon to avoid this step. */
				zval[m]+=(tz=exp(x)); inc=m*j;
				for(j=0;j<gmix->dimension;j++){
					mean[inc+j]+=(x=tz*feas->data[i][j]);
					dcov[inc+j]+=x*feas->data[i][j];
				}
			}
		}
	}
	workers_mutex_lock(t->mutex); /* Accumulate counts obtained to the mixture. */
	gmix->llh+=llh;
	for(m=0;m<gmix->num;m++){
		gmix->mix[m]._z+=zval[m]; inc=m*j;
		gmix->mix[m]._cfreq+=tfreq[m];
		for(j=0;j<gmix->dimension;j++){
			gmix->mix[m]._mean[j]+=mean[inc+j];
			gmix->mix[m]._dcov[j]+=dcov[inc+j];
		}
	}
	workers_mutex_unlock(t->mutex);
	free(zval); free(mean); free(tfreq);
}

/* Perform one iteration of the EM algorithm with the data and the mixture indicated. */
decimal gmm_EMtrain(data *feas,gmm *gmix,workers *pool){
	number m,i,inc,n=workers_number(pool); decimal tz,x;
	workers_mutex *mutex=workers_mutex_create();
	trainer *t=(trainer*)calloc(n,sizeof(trainer));
	/* Calculate expected value and accumulate the counts (E Step). */
	gmm_init_classifier(gmix);
	for(m=0;m<gmix->num;m++)
		gmix->mix[m]._cfreq=0;
	inc=feas->samples/n;
	for(i=0;i<n;i++){ /* Set and launch the parallel training. */
		t[i].feas=feas,t[i].gmix=gmix,t[i].mutex=mutex,t[i].ini=i*inc;
		t[i].end=(i==n-1)?(feas->samples):((i+1)*inc);
		workers_addtask(pool,thread_trainer,(void*)&t[i]);
	}
	workers_waitall(pool);
	workers_mutex_delete(mutex);
	/* Estimate the new parameters of the Gaussian Mixture (M Step). */
	for(m=0;m<gmix->num;m++){
		gmix->mix[m].prior=log((tz=gmix->mix[m]._z)/feas->samples);
		for(i=0;i<gmix->dimension;i++){
			gmix->mix[m].mean[i]=(x=gmix->mix[m]._mean[i]/tz);
			gmix->mix[m].dcov[i]=(gmix->mix[m]._dcov[i]/tz)-(x*x);
			if(gmix->mix[m].dcov[i]<gmix->mcov[i]) /* Smoothing covariances. */
				gmix->mix[m].dcov[i]=gmix->mcov[i];
		}
	}
	free(t);
	return gmix->llh/feas->samples;
}
