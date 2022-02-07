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

/* Allocate contiguous memory to create a new Gaussian Mixture. */
gmm *gmm_create(number n,number d){
	gmm *gmix=(gmm*)calloc(1,sizeof(gmm));
	gmix->mcov=(decimal*)calloc(gmix->dimension=d,sizeof(decimal));
	gmix->mix=(gauss*)calloc(gmix->num=n,sizeof(gauss));
	for(n=0;n<gmix->num;n++){
		gmix->mix[n].mean=(decimal*)calloc(gmix->dimension*4,sizeof(decimal));
		gmix->mix[n].dcov=gmix->mix[n].mean+gmix->dimension;
		gmix->mix[n]._mean=gmix->mix[n].dcov+gmix->dimension;
		gmix->mix[n]._dcov=gmix->mix[n]._mean+gmix->dimension;
	}
	return gmix;
}

/* Create and initialize the Mixture with maximum likelihood and disturb the means. */
gmm *gmm_initialize(data *feas,number nmix){
	gmm *gmix=gmm_create(nmix,feas->dimension);
	number i,j,k,b=feas->samples/gmix->num,bc=0;
	decimal x=1.0/gmix->num;
	/* Initialize the first Gaussian with maximum likelihood. */
	gmix->mix[0].prior=log(x);
	for(j=0;j<gmix->dimension;j++){
		gmix->mix[0].dcov[j]=x*feas->variance[j];
		gmix->mcov[j]=0.001*gmix->mix[0].dcov[j];
	}
	/* Disturb all the means creating C blocks of samples. */
	for(i=gmix->num-1;i>=0;i--,bc+=b){
		gmix->mix[i].prior=gmix->mix[0].prior;
		for(j=bc,x=0;j<bc+b;j++) /* Compute the mean of a group of samples. */
			for(k=0;k<gmix->dimension;k++)
				gmix->mix[i]._mean[k]+=feas->data[j][k];
		for(k=0;k<gmix->dimension;k++){ /* Disturbe the sample mean for each mixture. */
			gmix->mix[i].mean[k]=(feas->mean[k]*0.9)+(0.1*gmix->mix[i]._mean[k]/b);
			gmix->mix[i].dcov[k]=gmix->mix[0].dcov[k];
		}
	}
	return gmix;
}

/* Load the Gaussian Mixture from the file received as parameter. */
gmm *gmm_load(char *filename){
	number m,d; FILE *f=fopen(filename,"rb");
	if(!f)fprintf(stderr,"Error: Not %s model file found.\n",filename),exit(1);
	fread(&d,1,sizeof(number),f);
	fread(&m,1,sizeof(number),f);
	gmm *gmix=gmm_create(m,d);
	fread(gmix->mcov,sizeof(decimal)*gmix->dimension,1,f);
	for(m=0;m<gmix->num;m++){
		fread(&gmix->mix[m].prior,sizeof(decimal),1,f);
		fread(&gmix->mix[m].cgauss,sizeof(decimal),1,f);
		fread(gmix->mix[m].mean,gmix->dimension*sizeof(decimal)*2,1,f);
	}
	fclose(f);
	return gmix;
}

/* Save the Gaussian Mixture to the file received as parameter. */
void gmm_save(char *filename,gmm *gmix){
	number m; FILE *f=fopen(filename,"wb");
	if(!f)fprintf(stderr,"Error: Can not write to %s file.\n",filename),exit(1);
	fwrite(&gmix->dimension,1,sizeof(number),f);
	fwrite(&gmix->num,1,sizeof(number),f);
	fwrite(gmix->mcov,sizeof(decimal)*gmix->dimension,1,f);
	for(m=0;m<gmix->num;m++){
		fwrite(&gmix->mix[m].prior,sizeof(decimal),1,f);
		fwrite(&gmix->mix[m].cgauss,sizeof(decimal),1,f);
		fwrite(gmix->mix[m].mean,sizeof(decimal)*gmix->dimension*2,1,f);
	}
	fclose(f);
}

/* Save the trained model as a jSON log file (not usable). */
void gmm_save_log(char *filename,gmm *gmix){
	number m,j; FILE *f=fopen(filename,"w");
	if(!f)fprintf(stderr,"Error: Can not write to %s file.\n",filename),exit(1);
	fprintf(f,"{\n\t\"dimension\": %i,",gmix->dimension);
	fprintf(f,"\n\t\"classes\": %i,",gmix->num);
	fprintf(f,"\n\t\"minimum_dcov\": [ %.10f",gmix->mcov[0]);
	for(j=1;j<gmix->dimension;j++)
		fprintf(f,", %.10f",gmix->mcov[j]);
	fprintf(f," ],\n\t\"model\": [ ");
	for(m=0;m<gmix->num;m++){
		fprintf(f,"\n\t\t { \"class\": %i, \"lprior\": %.10f, ",m,gmix->mix[m].prior);
		fprintf(f,"\"means\": [ %.10f",gmix->mix[m].mean[0]);
		for(j=1;j<gmix->dimension;j++)
			fprintf(f,", %.10f",gmix->mix[m].mean[j]);
		fprintf(f," ], \"dcov\": [ %.10f",gmix->mix[m].dcov[0]);
		for(j=1;j<gmix->dimension;j++)
			fprintf(f,", %.10f",gmix->mix[m].dcov[j]);
		if(m==gmix->num-1)fprintf(f," ] }");
		else fprintf(f," ] },");
	}
	fprintf(f,"\n\t]\n}");
	fclose(f);
}

/* Free the allocated memory of the Gaussian Mixture. */
void gmm_delete(gmm *gmix){
	number i;
	for(i=0;i<gmix->num;i++)
		free(gmix->mix[i].mean);
	free(gmix->mcov);
	free(gmix->mix);
	free(gmix);
}
