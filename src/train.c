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

/* Show the trainer help message. */
void show_help(char *filename){
	fprintf(stderr,"Usage: %s <options>\n",filename);
	fprintf(stderr,"  Required:\n");
	fprintf(stderr,"    -d file.txt|file.gz   file that contains all the samples vectors\n");
	fprintf(stderr,"    -m file.gmm           file used to save the trained mixture model\n");
	fprintf(stderr,"  Recommended:\n");
	fprintf(stderr,"    -n 2-524228           optional number of components of the mixture\n");
	fprintf(stderr,"    -r file.json          optional file to save the log (not the model)\n");
	fprintf(stderr,"  Optional:\n");
	fprintf(stderr,"    -u 0.0-1.0            optional merge threshold based on similarity\n");
	fprintf(stderr,"    -s 0.0-1.0            optional stop criterion based on likelihood\n");
	fprintf(stderr,"    -i 1-9999             optional maximum number of EM iterations\n");
	fprintf(stderr,"    -t 1-256              optional maximum number of threads used\n");
	fprintf(stderr,"    -h                    optional argument that shows this message\n");
}

/* Show an error and leave the program. */
void show_error(const char *message){
	fprintf(stderr,"Error: %s.\n",message),exit(1);
}

/* Main execution of the trainer. */
int main(int argc,char *argv[]) {
	number i,o,x=0,nmix=-1,imax=99,t=sysconf(_SC_NPROCESSORS_ONLN);
	decimal last=INT_MIN,llh,sigma=0.1,m=-1.0; char *fnf=NULL,*fnm=NULL,*fnr=NULL;
	while((o=getopt(argc,argv,"r:u:t:i:d:m:n:s:h"))!=-1){
		switch(o){
			case 't': t=atoi(optarg);
				if(t>256||t<1)show_error("Number of threads must be on the 1-256 range");
				break;
			case 'n': nmix=atoi(optarg);
				if(nmix>524228||nmix<2)show_error("Number of components must be on the 2-32768 range");
				break;
			case 'i': imax=atoi(optarg);
				if(imax>9999||imax<1)show_error("Number of iterations must be on the 1-9999 range");
				break;
			case 'd': fnf=optarg,x++; break;
			case 'm': fnm=optarg,x++; break;
			case 'r': fnr=optarg; break;
			case 'u': m=atof(optarg);
				if(m>1.0||m<0.0)show_error("Merge threshold must be on the 0.0-1.0 range");
				break;
			case 's': sigma=atof(optarg);
				if(sigma>1.0||imax<0.0)show_error("Sigma criterion must be on the 0.0-1.0 range");
				break;
			case 'h': show_help(argv[0]),exit(1); break;
		}
	}
	if(x<2)show_help(argv[0]),exit(1); /* Test if exists all the needed arguments. */
	workers *pool=workers_create(t);
	data *feas=feas_load(fnf,pool); /* Load the features from the specified disc file. s*/
	if(nmix==-1){
		if(m<0)m=0.9;
		nmix=sqrt(feas->samples/2);
	}
	gmm *gmix=gmm_initialize(feas,nmix); /* Good GMM initialization using data.    */
	fprintf(stdout,"Components: %06i\n",gmix->num);
	for(o=1;o<=imax;o++){
		for(i=1;i<=imax;i++){
			llh=gmm_EMtrain(feas,gmix,pool); /* Compute one iteration of EM algorithm.   */
			fprintf(stdout,"%04i   %+7.2f%c   ln(likelihood)=(%.3f)\n",
				i,fabs(-100*(llh-last)/last),'%',llh); /* Show the EM results.  */
			if(last-llh>-sigma||isnan(last-llh))break; /* Break with sigma threshold. */
			last=llh;
		}
		x=gmix->num;
		if(m>=0){
			gmix=gmm_merge(gmix,feas,m,pool);
			fprintf(stdout,"Components: %06i\n",gmix->num);
		}
		if(x==gmix->num)break;
		last=INT_MIN;
	}
	workers_finish(pool);
	feas_delete(feas);
	if(fnr!=NULL)gmm_save_log(fnr,gmix);
	gmm_init_classifier(gmix); /* Pre-compute the non-data dependant part of classifier. */
	gmm_save(fnm,gmix); /* Save the model with the pre-computed part for fast classify.  */
	gmm_delete(gmix);
	return 0;
}
