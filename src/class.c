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

/* Show the classifier help message. */
void show_help(char *filename){
	fprintf(stderr,"Usage: %s <options>\n",filename);
	fprintf(stderr,"  Required:\n");
	fprintf(stderr,"    -d file.txt|file.gz   file that contains the samples vectors\n");
	fprintf(stderr,"    -m file.gmm           file of the trained model used to classify\n");
	fprintf(stderr,"  Recommended:\n");
	fprintf(stderr,"    -w file.gmm           optional world model used to smooth\n");
	fprintf(stderr,"    -r file.json          optional file to save the log (slower)\n");
	fprintf(stderr,"  Optional:\n");
	fprintf(stderr,"    -t 1-256              optional maximum number of threads used\n");
	fprintf(stderr,"    -h                    optional argument that shows this message\n");
}

/* Show an error and leave the program. */
void show_error(const char *message){
	fprintf(stderr,"Error: %s.\n",message),exit(1);
}

/* Main execution of the classifier. */
int main(int argc,char *argv[]){
	number i,o,x=0,t=sysconf(_SC_NPROCESSORS_ONLN); decimal result;
	char *fnr=NULL,*fnf=NULL,*fnm=NULL,*fnw=NULL;
	while((o=getopt(argc,argv,"t:d:m:w:r:h"))!=-1){
		switch(o){
			case 't': t=atoi(optarg);
				if(t>256||t<1)show_error("Number of threads must be on the 1-256 range");
				break;
			case 'r': fnr=optarg; break;
			case 'd': fnf=optarg,x++; break;
			case 'm': fnm=optarg,x++; break;
			case 'w': fnw=optarg; break;
			case 'h': show_help(argv[0]),exit(1); break;
		}
	}
	if(x<2)show_help(argv[0]),exit(1); /* Test if exists needed arguments. */
	workers *pool=workers_create(t);
	data *feas=feas_load(fnf,pool); /* Load the data from the specified file. */
	gmm *gmix=NULL,*gworld=NULL;
	if(fnw!=NULL)gworld=gmm_load(fnw); /* Load world model if is defined.  */
	gmix=gmm_load(fnm);
	if(fnr!=NULL)result=gmm_classify(fnr,feas,gmix,gworld,pool);
	else result=gmm_simple_classify(feas,gmix,gworld,pool);
	fprintf(stdout,"Score: %.10f\n",result);
	workers_finish(pool);
	gmm_delete(gmix);
	if(gworld!=NULL)gmm_delete(gworld);
	feas_delete(feas);
	return 0;
}
