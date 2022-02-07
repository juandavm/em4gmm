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
	decimal result;       /* Variable to store the result found. */
	data *feas;           /* Shared pointer to loaded samples.   */
	gmm *gmix;            /* Shared pointer to gaussian mixture. */
	gmm *gworld;          /* Shared pointer to gaussian mixture. */
	number ini, end;      /* Initial and final sample processed. */
	/* Data used only on detailed classifier to store the log.   */
	workers_mutex *mutex; /* Common mutex to lock shared data.   */
	FILE *f;              /* The jSON text file to save the log. */
	number *flag;         /* A flag for the first line on file.  */
}classifier;

/* Initialize the classifier by calculating the non-data dependant part. */
void gmm_init_classifier(gmm *gmix){
	decimal cache=gmix->dimension*log(2*NUM_PI);
	number m,j; gmix->llh=0;
	for(m=0;m<gmix->num;m++){
		gmix->mix[m].cgauss=gmix->mix[m]._z=0;
		for(j=0;j<gmix->dimension;j++){
			gmix->mix[m].cgauss+=log(gmix->mix[m].dcov[j]);
			gmix->mix[m].dcov[j]=1/gmix->mix[m].dcov[j];
			gmix->mix[m]._mean[j]=gmix->mix[m]._dcov[j]=0; /* Caches to 0. */
		}
		gmix->mix[m].cgauss=(2*gmix->mix[m].prior-(gmix->mix[m].cgauss+cache));
	}
}

/* Parallel and fast implementation of the Gaussian Mixture classifier. */
void thread_simple_classifier(void *tdata){
	classifier *t=(classifier*)tdata;
	decimal x,max,prob; number i,m,j;
	if(t->gworld!=NULL){ /* If the world model is defined, use it. */
		for(i=t->ini;i<t->end;i++){
			max=-HUGE_VAL;
			for(m=0;m<t->gworld->num;m++){
				prob=t->gworld->mix[m].cgauss;
				if(prob<max)continue; /* Speed-up the classifier. */
				for(j=0;j<t->gworld->dimension;j++){
					x=t->feas->data[i][j]-t->gworld->mix[m].mean[j];
					prob-=(x*x)*t->gworld->mix[m].dcov[j];
				}
				if(max<prob)max=prob;
			}
			t->result-=max; /* Compute final probability. */
		}
	}
	for(i=t->ini;i<t->end;i++){
		max=-HUGE_VAL;
		for(m=0;m<t->gmix->num;m++){
			prob=t->gmix->mix[m].cgauss; /* The precalculated non-data dependant part. */
			if(prob<max)continue; /* Speed-up the classifier. */
			for(j=0;j<t->gmix->dimension;j++){
				x=t->feas->data[i][j]-t->gmix->mix[m].mean[j];
				prob-=(x*x)*t->gmix->mix[m].dcov[j];
			}
			if(max<prob)max=prob; /* Fast classifier using Viterbi aproximation. */
		}
		t->result+=max;
	}
}

/* Efficient Gaussian Mixture classifier using a Viterbi aproximation. */
decimal gmm_simple_classify(data *feas,gmm *gmix,gmm *gworld,workers *pool){
	number i,n=workers_number(pool),inc=feas->samples/n; decimal result=0;
	classifier *t=(classifier*)calloc(n,sizeof(classifier));
	for(i=0;i<n;i++){ /* Set and launch the parallel classify. */
		t[i].feas=feas,t[i].gmix=gmix,t[i].gworld=gworld,t[i].ini=i*inc;
		t[i].end=(i==n-1)?(feas->samples):((i+1)*inc);
		workers_addtask(pool,thread_simple_classifier,(void*)&t[i]);
	}
	workers_waitall(pool); /* Wait to the end of the parallel classify. */
	for(i=0;i<n;i++)
		result+=t[i].result;
	free(t);
	return (result*0.5)/feas->samples;
}

/* Parallel implementation of the Gaussian Mixture classifier that holds the data. */
void thread_classifier(void *tdata){
	classifier *t=(classifier*)tdata;
	decimal x,max1,max2,prob; number i,m,j,c,s;
	char *buf1=(char*)calloc(s=25*t->gmix->num,sizeof(char));
	char *buf2=(char*)calloc(s=25*t->gmix->num,sizeof(char)),*bufi;
	for(i=t->ini;i<t->end;i++){
		if(t->gworld!=NULL){ /* If the world model is defined, use it. */
			max2=-HUGE_VAL;
			for(m=0;m<t->gworld->num;m++){
				prob=t->gworld->mix[m].cgauss;
				if(prob<max2)continue; /* Speed-up the classifier. */
				for(j=0;j<t->gworld->dimension;j++){
					x=t->feas->data[i][j]-t->gworld->mix[m].mean[j];
					prob-=(x*x)*t->gworld->mix[m].dcov[j];
				}
				if(max2<prob)max2=prob;
			}
		}else max2=0;
		/* Separated calculation of the first component to speed-up. */
		max1=t->gmix->mix[0].cgauss,c=0;
		for(j=0;j<t->gmix->dimension;j++){
			x=t->feas->data[i][j]-t->gmix->mix[0].mean[j];
			max1-=(x*x)*t->gmix->mix[0].dcov[j];
		}
		snprintf(buf1,s,"%.10f",(max1-max2)*0.5);
		for(m=1;m<t->gmix->num;m++){
			prob=t->gmix->mix[m].cgauss; /* The precalculated non-data dependant part. */
			for(j=0;j<t->gmix->dimension;j++){
				x=t->feas->data[i][j]-t->gmix->mix[m].mean[j];
				prob-=(x*x)*t->gmix->mix[m].dcov[j];
			}
			if(max1<prob)max1=prob,c=m; /* Fast classifier using Viterbi aproximation. */
			bufi=buf1; buf1=buf2; buf2=bufi;
			snprintf(buf1,s,"%s, %.10f",buf2,(prob-max2)*0.5);
		}
		t->result+=(max1-max2)*0.5;
		workers_mutex_lock(t->mutex); /* Write the classifier log on the jSON file. */
		t->gmix->mix[c]._cfreq++;
		if(t->flag[0]==0){
			fprintf(t->f,"\n\t\t{ \"sample\": %i, \"lprob\": [ %s ], \"class\": %i }",i,buf1,c);
			t->flag[0]=1;
		}else fprintf(t->f,",\n\t\t{ \"sample\": %i, \"lprob\": [ %s ], \"class\": %i }",i,buf1,c);
		workers_mutex_unlock(t->mutex);
	}
	free(buf1); free(buf2);
}

/* Detailed Gaussian Mixture classifier using a Viterbi aproximation. */
decimal gmm_classify(char *filename,data *feas,gmm *gmix,gmm *gworld,workers *pool){
	number i,n=workers_number(pool),inc=feas->samples/n; decimal result=0;
	workers_mutex *mutex=workers_mutex_create();
	classifier *t=(classifier*)calloc(n,sizeof(classifier));
	number *flag=(number*)calloc(1,sizeof(number));
	FILE *f=fopen(filename,"w");
	if(!f)fprintf(stderr,"Error: Can not write to %s file.\n",filename),exit(1);
	fprintf(f,"{\n\t\"samples\": %i,\n\t\"classes\": %i,",feas->samples,gmix->num);
	fprintf(f,"\n\t\"samples_results\": [ ");
	for(i=0;i<gmix->num;i++)
		gmix->mix[i]._cfreq=0;
	for(i=0;i<n;i++){ /* Set and launch the parallel classify. */
		t[i].feas=feas,t[i].gmix=gmix,t[i].gworld=gworld,t[i].ini=i*inc,t[i].mutex=mutex;
		t[i].end=(i==n-1)?(feas->samples):((i+1)*inc),t[i].f=f,t[i].flag=flag;
		workers_addtask(pool,thread_classifier,(void*)&t[i]);
	}
	workers_waitall(pool); /* Wait to the end of the parallel classify. */
	for(i=0;i<n;i++)
		result+=t[i].result;
	workers_mutex_delete(mutex);
	fprintf(f,"\n\t],\n\t\"mixture_occupation\": [ %i",gmix->mix[0]._cfreq);
	for(i=1;i<gmix->num;i++)
		fprintf(f,", %i",gmix->mix[i]._cfreq);
	fprintf(f," ],\n\t\"global_score\": %.10f\n}",result);
	fclose(f); free(t); free(flag);
	return result/feas->samples;
}
