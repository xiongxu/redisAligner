// gcc -O3 -Wall -W -Wstrict-prototypes -Wwrite-strings -g -ggdb -lgomp ref2aof.c -o ref2aof -I./hiredis -I./samtools-0.1.19 -L./hiredis -L./samtools-0.1.19 -lhiredis -lpthread -lz -lbam
// need 221G memory & over 12h
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <inttypes.h>
#include "IO_stream.h"
#include "faidx.h"
#include "dict.h"
#include "sds.h"

typedef struct BedGraph_t{
	char *chr;
	uint32_t start;
	uint32_t end;
	double score;
} BedGraph;

struct globalArgs_t {
	const char *infile;	/* -i option */
	const char *outfile;	/* -o option */
	const char *ref;	/* -r option */
} globalArgs;

void load_file(gzFile fp,const char *reference,const char *outfile);
void display_usage(char * argv[]);
static inline int readNextBedGraph(gzFile fz,char *buf,BedGraph *line);
static inline sds generateRedisCmd(sds *sdsArray);
static inline int str2int(const char *str);
static inline char *strupr(char *s);
static inline unsigned int myhashFunction(const void *key);
static inline void *mykeyDup(void *privdata, const void *key);
static inline int mykeyCompare(void *privdata, const void *key1, const void *key2);
static inline void mykeyDestructor(void *privdata, void *key);
static inline void myvalDestructor(void *privdata, void *obj);

uint8_t nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0/*'A'*/, 4, 1/*'C'*/,  4, 4, 4, 2/*'G'*/,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3/*'T'*/, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0/*'a'*/, 4, 1/*'c'*/,  4, 4, 4, 2/*'g'*/,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3/*'t'*/, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// static uint8_t groups[256] = {
//   1,7,7,0,2,10,6,8,3,11,8,5,10,8,6,6,
//   5,0,4,1,2,7,8,11,4,8,3,8,1,7,11,3,
//   4,3,4,4,7,6,6,7,1,5,1,6,0,5,9,8,
//   0,7,5,5,3,2,10,9,10,9,10,11,11,3,7,9,
//   7,3,5,11,3,2,6,8,8,11,10,7,3,8,6,7,
//   5,4,6,3,6,0,1,0,2,10,0,10,6,9,5,3,
//   1,5,2,6,10,6,9,5,3,2,2,1,10,9,6,11,
//   0,9,2,1,11,0,2,6,4,1,10,10,6,6,11,7,
//   9,11,0,9,8,11,6,6,7,5,2,0,6,6,8,3,
//   6,10,10,4,4,9,9,5,0,9,2,1,11,2,1,11,
//   1,7,3,4,8,5,7,1,11,11,9,10,8,9,2,5,
//   5,2,2,2,2,10,4,3,9,5,6,3,9,10,9,2,
//   2,11,6,7,5,9,2,1,1,5,1,8,10,9,10,4,
//   3,10,7,10,2,10,11,2,9,8,2,5,4,11,10,8,
//   1,4,9,2,11,0,8,3,2,0,8,0,3,5,3,4,
//   0,5,8,1,5,0,8,9,4,6,2,8,11,8,0,4
// };

static uint8_t groups[256] = {
	1,7,7,0,2,10,6,8,3,11,8,5,10,8,6,6,
	5,0,4,1,2,7,8,11,6,8,3,8,1,7,11,3,
	4,3,4,4,7,6,6,7,1,5,1,6,0,5,9,8,
	0,7,5,5,3,2,10,9,10,9,10,11,11,3,7,9,
	7,3,5,11,3,2,6,8,8,11,10,7,3,8,6,7,
	5,4,6,3,6,0,8,0,2,10,0,10,6,9,5,3,
	1,5,2,6,10,6,9,5,3,2,2,1,10,9,6,11,
	0,9,2,1,11,0,2,6,4,1,10,10,6,6,11,7,
	9,11,0,9,8,11,6,6,7,5,2,0,6,6,8,3,
	6,10,10,4,4,9,9,5,0,9,2,1,11,2,1,11,
	1,7,3,4,8,5,7,1,11,11,9,10,8,9,2,5,
	5,2,2,2,2,10,6,3,9,5,6,3,9,10,9,2,
	2,11,6,7,5,9,2,1,1,5,1,8,10,9,10,4,
	3,10,7,10,2,10,11,2,9,8,2,5,4,11,10,8,
	1,4,9,2,11,0,8,3,2,0,8,0,3,5,3,4,
	0,5,8,1,5,0,8,9,4,6,2,8,11,8,0,4
};

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c) InfiniGenomics 2017\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for load & pack the uniquely mapped region into redis \n" \
"Usage: %s [-i Infile] -r [reference]\n" \
"Example1:\n  %s -i ../redisApp/wgEncodeCrgMapabilityAlign36mer_1.BedGraph -o hg19 -r ../redisApp/hg19.fasta  \n" \
"\n" \
"   [-i Infile]    = Infile.default is stdin                            [required]\n" \
"   [-r ref]       = reference.                                         [required]\n" \
"   [-h] This helpful help screen.                                      [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

#define InsertHash(hash,key,key2,value,mydictType) \
do{																\
	dictEntry *entry=dictFind(hash,key),*entry2=NULL;			\
	dict *hashChild=NULL;										\
	if (entry==NULL){											\
		hashChild=dictCreate(mydictType,NULL);					\
		dictAdd(hashChild,key2,value);							\
		dictAdd(hash,key,hashChild);							\
	}else{														\
		hashChild=(dict *)entry->val;							\
		entry2=dictFind(hashChild,key2);						\
		if (entry2==NULL){										\
			dictAdd(hashChild,key2,value);						\
		}else{													\
			fprintf(stderr,"replace\n");						\
			dictReplace(hashChild,key2,value);					\
		}														\
	}															\
}while(0)

#define initdictType(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=NULL;						\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=myvalDestructor;	\
}while (0)

#define initdictType2(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=mykeyDup;					\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=mykeyDestructor;	\
}while (0)

unsigned int myhashFunction(const void *key){
	return dictGenHashFunction((unsigned char *)key,sdslen((sds) key));
}

void *mykeyDup(void *privdata, const void *key){
	return (void *)sdsdup((sds) key);
}

int mykeyCompare(void *privdata, const void *key1, const void *key2){
	return !sdscmp((sds) key1,(sds) key2);
}

void mykeyDestructor(void *privdata, void *key){
	sdsfree((sds) key);
}

void myvalDestructor(void *privdata, void *obj){
	dict *dictObj=(dict *)obj;
	dictRelease(dictObj);
}

int readNextBedGraph(gzFile fz,char *buf,BedGraph *line){
	buf=gzgets(fz,buf,512);
	if (gzeof(fz)) return 0;
	// BedGraph *line=(BedGraph *)malloc(sizeof(BedGraph));
	// line->chr=(char *)calloc(64,sizeof(char));
	sscanf(buf, "%s\t%u\t%u\t%lf\n",line->chr,&line->start,&line->end,&line->score);
//	line->chr=realloc(line->chr,(strlen(line->chr)+1)*sizeof(char));
//	fprintf(stderr,"%s\t%u\t%u\t%lf\n",line->chr,line->start,line->end,line->score);
	return 1;
}

int str2int(const char *str){
	int temp = 0;
	const char *ptr = str;
	if (*ptr == '-' || *ptr == '+') ptr++;
	while(*ptr){
		if (*ptr == 'X' || *ptr == 'x') return 23;
		if (*ptr == 'Y' || *ptr == 'y') return 24;
		if (*ptr == 'M' || *ptr == 'm') return 25;
		if ((*ptr < '0') || (*ptr > '9')) {ptr++;continue;}
		temp = temp * 10 + (*ptr - '0');
		ptr++;
	}
	if (*ptr == '-') temp = -temp;
	return temp;
}

char *strupr(char *s){
	char *p=s;
	while ((*p++=toupper(*p))!=0); //toupper was included in ctype.h
	return s;
}

sds generateRedisCmd(sds *sdsArray){
	sds s=sdsnewlen("*4\r\n$4\r\nHSET\r\n",14);
	int i,Length=0;
	for (i=0;i<3;++i) {
		Length=sdslen(sdsArray[i]);
		s=sdscatprintf(s,"$%d\r\n",Length);
		s=sdscatlen(s,sdsArray[i],Length);
		s=sdscatlen(s,"\r\n",2);
	}
	return s;
}

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

void load_file(gzFile fp,const char *reference,const char *outfile) {
	long long begin,n=0;
	begin=usec();
	faidx_t *ref_fai=fai_load(reference);
	BedGraph *line=(BedGraph *)malloc(sizeof(BedGraph));
	line->chr=(char *)calloc(64,sizeof(char));
	int seq_len=0,i=0;
	uint32_t p=0;
	char *region=(char *)calloc(64,sizeof(char)),*buf = (char *) calloc(512,sizeof(char)),*seq=NULL;
	sds *sdsArray=(sds *)malloc(3*sizeof(sds));
	sdsArray[0]=sdsnewlen(NULL,3);
	sdsArray[1]=sdsnewlen(NULL,6);
	sdsArray[2]=sdsnewlen(NULL,5);
	
	dictType *mydictType=(dictType *)malloc(sizeof(dictType));
	initdictType(mydictType);
	dictType *mydictType2=(dictType *)malloc(sizeof(dictType));
	initdictType2(mydictType2);

	dict *ht=dictCreate(mydictType,NULL);
	uint64_t stat[256]={0},groupStat[12]={0},count=0;
	FILE **aofs=(FILE **)malloc(12*sizeof(FILE *));
	for(i=0;i<12;++i){
		sprintf(region,"%d.aof",i);
		aofs[i]=fcreat_outfile(outfile,region);
		fprintf(aofs[i],"*2\r\n$6\r\nSELECT\r\n$1\r\n0\r\n");
	}
	while (readNextBedGraph(fp,buf,line)) {
		if (line->score <1 || !strcmp(line->chr,"chrM")) continue;
		for (p=line->start+1;p<=line->end ;p++ ){
			sprintf(region,"%s:%d-%d",line->chr,p,p+35);
			seq=fai_fetch(ref_fai,region,&seq_len);
			memset(sdsArray[0],0,3);
			memset(sdsArray[1],0,6);
			memset(sdsArray[2],0,5);
			for(i=0;i<12;++i){
				_set_pac((uint8_t *)sdsArray[0],i,nst_nt4_table[(uint8_t)seq[i]]);
			}
			groupStat[groups[(uint8_t)sdsArray[0][0]]]++;
			stat[(uint8_t)sdsArray[0][0]]++;
			n++;
			for(;i<seq_len;++i){
				_set_pac((uint8_t *)sdsArray[1],i-12,nst_nt4_table[(uint8_t)seq[i]]);
			}
			uint8_t GCcount=0;
			for(i=0;i<seq_len;++i){
				if (nst_nt4_table[(uint8_t)seq[i]]==1 || nst_nt4_table[(uint8_t)seq[i]]==2) GCcount++;
			}
			sdsArray[2][0] = (uint8_t)str2int(line->chr) | ((GCcount&0x70)<<1);
			sdsArray[2][1] = ((GCcount&0xf)<<4) | (p>>24);
			for (i=0;i<=2 ;++i){
				sdsArray[2][4-i] = p>>(i*8) & 0xff ;
			}
			InsertHash(ht,sdsArray[0],sdsArray[1],sdsArray[2],mydictType2);
			free(seq);
		}
	}
	free(line->chr);
	free(line);
	for(i=0;i<3;++i){
		sdsfree(sdsArray[i]);
	}
	free(sdsArray);
	dictIterator *iter=dictGetIterator(ht),*iter2=NULL;
	dictEntry *entry=dictNext(iter),*entry2=NULL;
	for (; entry; entry=dictNext(iter)) {
		FILE *aof=aofs[groups[(uint8_t)((sds)entry->key)[0]]];
		dict *hashChild=(dict *)entry->val;
		sds CMD=sdsnewlen("*",1);
		int key1Len=sdslen((sds)entry->key);
		CMD=sdscatprintf(CMD,"%lu\r\n$5\r\nHMSET\r\n$%d\r\n",2*(1+hashChild->used),key1Len);
		CMD=sdscatlen(CMD,(sds)entry->key,key1Len);
		CMD=sdscatlen(CMD,"\r\n",2);
		iter2=dictGetIterator(hashChild);
		for (entry2=dictNext(iter2); entry2; entry2=dictNext(iter2)) {
			int key2Len=sdslen((sds)entry2->key);
			CMD=sdscatprintf(CMD,"$%d\r\n",key2Len);
			CMD=sdscatlen(CMD,(sds)entry2->key,key2Len);
			int valueLen=sdslen((sds)entry2->val);
			CMD=sdscatprintf(CMD,"\r\n$%d\r\n",valueLen);
			CMD=sdscatlen(CMD,(sds)entry2->val,valueLen);
			CMD=sdscatlen(CMD,"\r\n",2);
		}
		fwrite(CMD,1,sdslen(CMD),aof);
		sdsfree(CMD);
	}
	dictReleaseIterator(iter);
	dictRelease(ht);
	free(mydictType);
	free(mydictType2);
	for(i=0;i<256;++i){
		fprintf(stderr,"%d\t%hhu\t%" PRIu64 "\n",i,groups[i],stat[i]);
		count+=stat[i];
	}
	fprintf(stderr,"\nkey count: %" PRIu64 "\n\n",count);
	for(i=0;i<12;++i){
		fclose(aofs[i]);
		fprintf(stderr,"%d\t%" PRIu64 "\n",i,groupStat[i]);
	}
	fprintf(stderr,"Finished %lld in %.3f s\n",n,(double)(usec()-begin)/CLOCKS_PER_SEC);
	free(buf);
	free(region);
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.ref="-";
	globalArgs.outfile="-";
	const char *optString = "i:r:o:h?";
	if (argc<2){
		display_usage(argv);
		exit(1);
	}
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'i':
				globalArgs.infile = optarg;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'r':
				globalArgs.ref = optarg;
				break;
			case '?':	/* fall-through is intentional */
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				/* You won't actually get here. */
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	gzFile in=open_input_stream(globalArgs.infile);
	load_file(in,globalArgs.ref,globalArgs.outfile);
	gzclose(in);
	return 0;
}
