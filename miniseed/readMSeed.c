/*-----------------------------------------------------------------------
 Test code to read a miniSEED file
------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "libmseed.h"

MSTrace*
readTraceMSeed(char* filename, char* sdt, char* edt, int* nsamp) 
{
	MSRecord* msr = NULL;
	MSFileParam* msfp = NULL;
	MSTrace* mst;
	int retcode,rettrace,nrec;
	flag verbose = 0;
	int jday,year,month,mday,hour,min,sec,nsec;
	hptime_t tanf,tend;
	char isotimestr[27];

	sscanf(sdt,"%4d%2d%2d_%2d%2d%2d_%9d",&year,&month,&mday,&hour,&min,&sec,&nsec);  /* handle start time */
	ms_md2doy(year,month,mday,&jday);
	tanf = ms_time2hptime(year,jday,hour,min,sec,nsec/1000);	

	sscanf(edt,"%4d%2d%2d_%2d%2d%2d_%9d",&year,&month,&mday,&hour,&min,&sec,&nsec); /* handle end time */
	ms_md2doy(year,month,mday,&jday);
	tend = ms_time2hptime(year,jday,hour,min,sec,nsec/1000);

	mst = malloc(sizeof(MSTrace));
	memset (mst, 0, sizeof (MSTrace));
	mst->numsamples = 0;
	mst->samplecnt = 0;
	nrec = 0;
	while ((retcode = ms_readmsr_r(&msfp, &msr,filename,0,NULL,NULL,1,1,verbose)) == MS_NOERROR)
	{
/*		printf("Record start time: %s\n",ms_hptime2isotimestr(msr->starttime,isotimestr,1));
		printf("Record end   time: %s\n",ms_hptime2isotimestr(msr_endtime(msr),isotimestr,1));
		printf("Tanf             : %s\n",ms_hptime2isotimestr(tanf,isotimestr,1));
		printf("Tend             : %s\n",ms_hptime2isotimestr(tend,isotimestr,1)); */
		if ((msr->starttime < tanf) && (msr_endtime(msr) < tanf)) continue;     /* skip record */
		if (msr->starttime > tend) break;                                       /* exit loop */
		nrec++;                                                                 /* record fits time window */
		if (nrec == 1) {                                                        /* fill part of mst-structure */
			ms_strncpclean(mst->network,msr->network,11);
			ms_strncpclean(mst->station,msr->station,11);
			ms_strncpclean(mst->location,msr->location,11);
			ms_strncpclean(mst->channel,msr->channel,11);
			mst->dataquality = msr->dataquality;
			mst->starttime = msr->starttime;
			mst->samprate = msr->samprate;
			mst->sampletype = msr->sampletype;
		}
/*		printf("Record start time: %s\n",ms_hptime2isotimestr(msr->starttime,isotimestr,1)); */
		rettrace = mst_addmsr(mst,msr,1);                 /* append data record to trace structure */
		if (rettrace == -1) 
		{
			retcode = MS_GENERROR;
			break;
		}
/*		printf("Trace end time:    %s\n",ms_hptime2isotimestr(mst->endtime,isotimestr,1)); */
	}
	retcode = ms_readmsr_r(&msfp,&msr,NULL,0,NULL,NULL,1,0,verbose);
	if (nrec == 0) {
		*nsamp = 0;
		return NULL;
	}
	if ((retcode != MS_NOERROR) && (retcode != MS_ENDOFFILE)) {    /* error when reading file */
		ms_log (2, "Error reading input file %s: %s\n", filename, ms_errorstr(retcode));
		*nsamp = 0;
		return NULL;
	}
	*nsamp = mst->numsamples;
	return mst;
}

MSTraceGroup*
readTraceGroupMSeed(char* filename, int* numtraces) 
{
	MSRecord* msr = NULL;
	MSFileParam* msfp = NULL;
	MSTrace* mst = NULL;
	MSTraceGroup *mstg = NULL;
	int retcode,rettrace,nrec;
	flag verbose = 0;

	mstg = malloc(sizeof(MSTraceGroup));
	memset (mstg, 0, sizeof (MSTraceGroup));
	while ((retcode = ms_readmsr_r(&msfp, &msr,filename,0,NULL,NULL,1,1,verbose)) == MS_NOERROR)
	{
/*		msr_print(msr, verbose-1);       */
/*
 * Reset MSRecord data samples to avoid accumulating all samples in the MSTrace 
 */
		mst = mst_addmsrtogroup (mstg, msr, 0, -1.0, -1.0);
		msr->datasamples = 0;
		msr->numsamples = 0;
	}
	retcode = ms_readmsr_r(&msfp,&msr,NULL,0,NULL,NULL,1,0,verbose);
	if ((retcode != MS_NOERROR) && (retcode != MS_ENDOFFILE)) {    /* error when reading file */
		ms_log (2, "Error reading input file %s: %s\n", filename, ms_errorstr(retcode));
		return NULL;
	}
	retcode = mst_groupsort(mstg,0);
	if (retcode < 0) {
		ms_log (2, "Sorting unsuccessful: %s\n", filename, ms_errorstr(retcode));
		return NULL;
	}
	*numtraces = mstg->numtraces;
	return mstg;
}

MSTrace*
getFirstTraceMSeed(MSTraceGroup* mstg, int* nsamp) 
{
    	*nsamp = (mstg->traces)->numsamples;
	return mstg->traces;
}

MSTrace*
getNextTraceMSeed(MSTrace* mstcur, int* nsamp) 
{
    	MSTrace* next;
	next = mstcur->next;
	if (next) {
	    	*nsamp = next->numsamples;
	} else {
	    	*nsamp = 0;
	}
	return next;
}

int
unpackTraceMSeed(MSTrace* mst,int* n,double* dt,char* sdate,int* tfs,int* tns,
	         char* network,char* station,char* channel,char* location,float y[],int* mstfree)
{
	int* idata;
	float* fdata;
	double* ddata;
	int year,month,mday,hour,min,sec,usec,i;
	char isotimestr[27];
	char z[1];

/*	printf("Station %s, Channel %s, Network %s, Location %s\n",mst->station,mst->channel,mst->network,mst->location);
	printf("Sample type is %c\n",mst->sampletype);
	printf("True start time: %s\n",ms_hptime2isotimestr(mst->starttime,isotimestr,1));
	printf("True end time:   %s\n",ms_hptime2isotimestr(mst->endtime,isotimestr,1));
	printf("Samples %d\n",mst->numsamples); */
	*n = mst->numsamples;
	*dt = 1./mst->samprate;
	ms_hptime2isotimestr(mst->starttime,isotimestr,1);     /* extract start time in detail */
	sscanf(isotimestr,"%4d%c%2d%c%2d%c%2d%c%2d%c%2d%c%6d",&year,z,&month,z,&mday,z,&hour,z,&min,z,&sec,z,&usec);
	sprintf(sdate,"%4.4d%2.2d%2.2d",year,month,mday);
	*tfs = (hour*60+min)*60+sec;                           /* start time to full seconds after midnight */
	*tns = usec*1000;                                      /* fractional part in nanoseconds */
	ms_strncpclean(network,mst->network,11);               /* set station etc */
	ms_strncpclean(station,mst->station,11);
	ms_strncpclean(channel,mst->channel,11);
	ms_strncpclean(location,mst->location,11);
	if (mst->sampletype == 'i') {                          /* convert to float from sample type */
		idata = (int* ) mst->datasamples;
		for (i=0; i < *n; i++) {y[i] = (float) idata[i];}
	} else if(mst->sampletype == 'f') {
		fdata = (float* ) mst->datasamples;
		for (i=0; i < *n; i++) {y[i] = fdata[i];}
	} else if(mst->sampletype == 'd') {
		ddata = (double* ) mst->datasamples;
		for (i=0; i < *n; i++) {y[i] = (float) ddata[i];}
	} else {
		ms_log (2, "Non-usable datatype %c: %s\n", mst->sampletype, ms_errorstr(MS_GENERROR));
		return -1;
	}	
	if (*mstfree) mst_free(&mst);                                        /* free trace */
	return 0;
}
MSTrace*
readFullTraceMSeed(char* filename, int* nsamp) 
{
	MSRecord* msr = NULL;
	MSFileParam* msfp = NULL;
	MSTrace* mst;
	int retcode,rettrace,nrec;
	flag verbose = 0;
	hptime_t tanf,tend;

	mst = malloc(sizeof(MSTrace));
	memset (mst, 0, sizeof (MSTrace));
	mst->numsamples = 0;
	mst->samplecnt = 0;
	nrec = 0;
	while ((retcode = ms_readmsr_r(&msfp, &msr,filename,0,NULL,NULL,1,1,verbose)) == MS_NOERROR)
	{
		nrec++;                                                                 /* record fits time window */
		if (nrec == 1) {                                                        /* fill part of mst-structure */
			ms_strncpclean(mst->network,msr->network,11);
			ms_strncpclean(mst->station,msr->station,11);
			ms_strncpclean(mst->location,msr->location,11);
			ms_strncpclean(mst->channel,msr->channel,11);
			mst->dataquality = msr->dataquality;
			mst->starttime = msr->starttime;
			mst->samprate = msr->samprate;
			mst->sampletype = msr->sampletype;
		}
		rettrace = mst_addmsr(mst,msr,1);                 /* append data record to trace structure */
		if (rettrace == -1) 
		{
			retcode = MS_GENERROR;
			break;
		}
	}
	retcode = ms_readmsr_r(&msfp,&msr,NULL,0,NULL,NULL,1,0,verbose);
	if (nrec == 0) {
		*nsamp = 0;
		return NULL;
	}
	if ((retcode != MS_NOERROR) && (retcode != MS_ENDOFFILE)) {    /* error when reading file */
		ms_log (2, "Error reading input file %s: %s\n", filename, ms_errorstr(retcode));
		*nsamp = 0;
		return NULL;
	}
	*nsamp = mst->numsamples;
	return mst;
}
/*------------------------------------------------------------------------
 *  pack a MS record into a miniSeed file
 -------------------------------------------------------------------------*/
static void 
record_handler (char *record, int reclen, void *srcname) 
{
	if ( fwrite(record, 1, reclen, srcname) == 0 ) {
		ms_log (2, "Error writing to output file");
	}
}

int 
writeTraceMSeed(char* filename, int* idata, int* nsamp, double* samprate, char* timestr,
	          char* network, char* station, char* channel, char* location)
{	
	MSRecord* msr;
	FILE* fpms;
	int psamples;
	int precords, i;
/**/
	printf("writeTraceMSeed: filename = %s, nsamp = %d, samprate = %f, timestr = %s\n",filename,*nsamp,*samprate,timestr);

/*  fill MSTrace structure  */

	printf("Some data: %d,%d,%d,%d,%d\n",idata[0],idata[1],idata[2],idata[3],idata[4]);
	msr = msr_init(NULL);
	strcpy(msr->network, network);
	strcpy(msr->station, station);
	strcpy(msr->channel, channel);
	strcpy(msr->location, location);
	msr->dataquality = 'D';
	msr->starttime = ms_seedtimestr2hptime(timestr);
	msr->samprate = *samprate;
	msr->datasamples = (void *) idata;
	msr->numsamples = *nsamp;
	msr->sampletype = 'i';
	msr->reclen = 4096;
	msr->byteorder = 1;
	msr->encoding = DE_STEIM2;

/*  open output file and write record */

	fpms = fopen(filename, "w");
	precords = msr_pack (msr, &record_handler, fpms, &psamples, 1, 0);
	ms_log (0, "Packed %d samples into %d records\n",psamples, precords);

/* clean up */

	msr_free (&msr);
	fclose(fpms);
	return precords;
}
