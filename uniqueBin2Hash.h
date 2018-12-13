//gcc -g -O3 -Wall get_cv_outlier.c -o get_cv_outlier -I/home/xuxiong/work/c/seq2sds/hiredis-master -I/home/xuxiong/work/c -I/share/software/software/zlib_1.2.8_install/include -L/share/software/software/zlib_1.2.8_install/lib -L/home/xuxiong/work/c/seq2sds/hiredis-master -lz -lm -lhiredis
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "IO_stream.h"
#include "kvec.h"
#include "dict.h"

#define TABLESIZE 50000
#define FreeBedGraph(line) \
do{							\
	free((line)->chr);		\
	free((line));			\
}while (0)

typedef struct _BedGraph_ {
	char *chr;
	uint32_t start;
	uint32_t end;
	double score;
	uint32_t length;
	uint32_t sumLen;
	uint32_t bin;
} BedGraph;

typedef struct _bin_ {
	uint32_t bin;
	uint32_t start;
	uint32_t end;
} Bin;

typedef struct _array_ {
	Bin *bins;
	int alloc;
	int count;
} ArrayBin;

typedef kvec_t(Bin) vecBin;
typedef kvec_t(vecBin) vecVecBin;

static inline long long usec(void);
static inline BedGraph *readNextBedGraph(char *buf);
dictType *initDictType(void);
void outputSortedHash(dict *ht,char *outfile);
dictEntry** dump_dict(dict *ht);
static inline int compare_hashed_key(const void *a, const void *b);
static inline int str2int(const char *str);
dict *load_bin(const char *infile,dictType *mydictType,int BinLen);
void output_hash(dict *hashtbl,char *outfile);
static inline void overlap(int *j,uint32_t chrLen ,int windows,vecBin *indexedBins,Bin *bin,uint16_t windowSize);
void IndexBins(dict *ht,uint16_t unitLen);
//vecBin **indexUniqueContinousBin(vecBin *chromosomesContinousBin,int unitLen);
//void outputIndexedUniqueContinousBin(vecBin **indexedBins,vecBin *chromosomesContinousBin,int unitLen);
void outputIndexedUniqueContinousBin(vecVecBin *indexedBins,vecBin *chromosomesContinousBin,int unitLen);
vecVecBin *indexUniqueContinousBin(vecBin *chromosomesContinousBin,uint16_t unitLen);

vecBin *loadUnique2BinFile(const char *infile);

#define initdictType(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=NULL;						\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=myvalDestructor;	\
}while (0)

static inline void *mykeyDup(void *privdata, const void *key);
static inline int mykeyCompare(void *privdata, const void *key1, const void *key2);
static inline void mykeyDestructor(void *privdata, void *key);
static inline void myvalDestructor(void *privdata, void *obj);
static inline unsigned int myhashFunction(const void *key);

unsigned int myhashFunction(const void *key){
	return dictGenHashFunction((unsigned char *)key,strlen((char *)key));
}

void *mykeyDup(void *privdata, const void *key) {
	return (void *)strdup((char *) key);
}

int mykeyCompare(void *privdata, const void *key1, const void *key2) {
	size_t l1,l2,minlen;
	int cmp;
	l1 = strlen((char *)key1);
	l2 = strlen((char *)key2);
	minlen = (l1 < l2) ? l1 : l2;
	cmp = memcmp((char *)key1,(char *)key2,minlen);
	if (cmp == 0) return !(l1-l2);
	return !cmp;
}

void mykeyDestructor(void *privdata, void *key){
	if (key == NULL) return;
	free(key);
}

void myvalDestructor(void *privdata, void *obj){
	ArrayBin *myobj=(ArrayBin *)obj;
	free(myobj->bins);
	free(myobj);
}

BedGraph *readNextBedGraph(char *buf) {
	BedGraph *line = (BedGraph *)malloc(sizeof(BedGraph));
	line->chr=(char *)calloc(16,sizeof(char));
	sscanf(buf,"%s\t%u\t%u\t%lf\t%u\t%u\t%u[^\n]",line->chr,&(line->start),&(line->end),&(line->score),&(line->length),&(line->sumLen),&(line->bin));
	line->chr=realloc(line->chr,(strlen(line->chr)+1)*sizeof(char));
	return line;
}

dictType *initDictType(void){
	dictType *mydictType=(dictType *)malloc(sizeof(dictType));
	initdictType(mydictType);
	return mydictType;
}

dict *load_bin(const char *infile,dictType *mydictType,int BinLen){
	dict *ht=dictCreate(mydictType,NULL);
	gzFile fq=open_input_stream(infile);
	BedGraph *last_line=NULL;
	uint32_t i=0;
	char *buf=(char *)calloc(1024,sizeof(char));
	while (gzgets(fq,buf,1024*sizeof(char))!=NULL){
		BedGraph *line=readNextBedGraph(buf);
		dictEntry *entry=dictFind(ht,line->chr);
		if (entry==NULL){
			if (last_line){
				entry=dictFind(ht,last_line->chr);
				if (entry) {
					ArrayBin *current_array=(ArrayBin *)entry->val;
					(current_array->bins+current_array->count++)->end = last_line->end;
				}
			}
			ArrayBin *array=malloc(sizeof(ArrayBin));
			array->alloc=TABLESIZE;
			array->count=0;
			array->bins=(Bin *)calloc(array->alloc,sizeof(Bin));
			(array->bins+array->count)->start=line->start;
			(array->bins+array->count)->bin=line->bin;
			dictAdd(ht, line->chr,array);
		}else{
			ArrayBin *current_array=(ArrayBin *)entry->val;
			if (line->bin != last_line->bin) {
				if (line->bin-last_line->bin>1){
					for (i=1; i<=(last_line->sumLen % BinLen +line->length)/BinLen ;i++ ){
						uint32_t current_end=line->start+ i*BinLen - (last_line->sumLen % BinLen);
						(current_array->bins+current_array->count++)->end = current_end;
						//(current_array->bins+current_array->count)->start = current_end+1;
						(current_array->bins+current_array->count)->start = current_end;
						(current_array->bins+current_array->count)->bin = last_line->bin+i;
					}
				}else{
					uint32_t current_end=line->end-(line->sumLen % BinLen);
					(current_array->bins+current_array->count++)->end = current_end;
					//(current_array->bins+current_array->count)->start = current_end+1;
					(current_array->bins+current_array->count)->start = current_end;
					(current_array->bins+current_array->count)->bin = line->bin;
				}
			}
		}
		if (last_line) {
			FreeBedGraph(last_line);
		}
		last_line=line;
	}
	if (last_line){
		dictEntry *entry=dictFind(ht,last_line->chr);
		if (entry) {
			ArrayBin *current_array=(ArrayBin *)entry->val;
			(current_array->bins+current_array->count)->end = last_line->end;
		}
		FreeBedGraph(last_line);
	}
	free(buf);
	gzclose(fq);
	return ht;
}

int str2int(const char *str){
	int temp = 0;
	const char *ptr = str;
	if (*str == '-' || *str == '+') str++;
	while(*str){
		if (*str == 'X') {
			return 23;
		}else if (*str == 'Y') {
			return 24;
		}else if (*str == 'M') {
			return 25;
		}
		if ((*str < '0') || (*str > '9')) {str++;continue;}
		temp = temp * 10 + (*str++ - '0');
	}
	if (*ptr == '-') temp = -temp;
	return temp;
}

int compare_hashed_key(const void *a, const void *b) {
	return str2int((char *)(*(dictEntry* const *)a)->key)-str2int((char *)(*(dictEntry* const *)b)->key);
}

dictEntry** dump_dict(dict *ht){
	dictEntry** D = (dictEntry **)malloc(ht->used * sizeof(dictEntry *));
	dictIterator *iter=dictGetIterator(ht);
	dictEntry *entry=dictNext(iter);
	uint32_t j=0;
	for (; entry; entry=dictNext(iter)) {
		D[j++]=entry;
	}
	dictReleaseIterator(iter);
	qsort(D, ht->used, sizeof(dictEntry *), compare_hashed_key);
	return D;
}

void output_hash(dict *hashtbl,char *outfile){
	FILE *out1=fcreat_outfile(outfile,".bin2pos");
	dictIterator *iter=dictGetIterator(hashtbl);
	dictEntry *entry=dictNext(iter);
	int i;
	for (;entry;entry=dictNext(iter) ){
		ArrayBin *array=(ArrayBin *)(entry->val);
		for (i=0;i<array->count ;i++ ){
			fprintf(out1,"%s\t%u\t%u\t%u\n",(char *)entry->key,array->bins[i].start,array->bins[i].end,array->bins[i].bin);
		}
	}
	dictReleaseIterator(iter);
	fclose(out1);
}

void outputSortedHash(dict *ht,char *outfile){
	if (ht){
		FILE *out=fcreat_outfile(outfile,".bin2pos");
		dictEntry** D=dump_dict(ht);
		unsigned long chri=0;
		int i;
		for (chri=0;chri<ht->used;chri++) {
			ArrayBin *array=(ArrayBin *)(D[chri]->val);
			for (i=0;i<array->count ;i++ ){
				fprintf(out,"%s\t%u\t%u\t%u\n",(char *)D[chri]->key,array->bins[i].start,array->bins[i].end,array->bins[i].bin);
			}
		}
		free(D);
		fclose(out);
	}
}

vecBin *loadUnique2BinFile(const char *infile) {
	gzFile fq=open_input_stream(infile);
	char *buf=(char *)calloc(1024,sizeof(char)),currentChr[32]={0};
	vecBin *chromosomesContinousBin = (vecBin *)calloc(25,sizeof(vecBin));
	while (gzgets(fq,buf,1024*sizeof(char))!=NULL) {
		Bin *line = (Bin *)calloc(1,sizeof(Bin));
		sscanf(buf,"%s\t%u\t%u\t%u[^\n]",currentChr,&(line->start),&(line->end),&(line->bin));
		int chrIndex = str2int(currentChr)-1;
		if(chromosomesContinousBin[chrIndex].m){
			kv_push(Bin, chromosomesContinousBin[chrIndex], *line);
		}else{
			kv_init(chromosomesContinousBin[chrIndex]);
			kv_push(Bin, chromosomesContinousBin[chrIndex], *line);
		}
	}
	free(buf);
	gzclose(fq);
	fprintf(stderr,"Done loadUnique2BinFile\n");
	return chromosomesContinousBin;
}

vecVecBin *indexUniqueContinousBin(vecBin *chromosomesContinousBin,uint16_t unitLen) {
	int i=0;
	vecVecBin *indexedBins = (vecVecBin *)calloc(25,sizeof(vecVecBin));
	for(i=0;i<25;++i){
		uint32_t chrEnd = chromosomesContinousBin[i].a[chromosomesContinousBin[i].n-1].end;
		int window_count = chrEnd/unitLen+1;
		indexedBins[i].a = (vecBin *)calloc(window_count, sizeof(vecBin));
		indexedBins[i].m = indexedBins[i].n = window_count;
		int j=0; size_t k=0;
		for (k=0;k<chromosomesContinousBin[i].n ;++k ) {
			overlap(&j,chrEnd,window_count,indexedBins[i].a,chromosomesContinousBin[i].a+k,unitLen);
		}
	}
	fprintf(stderr,"Done indexUniqueContinousBin\n");
	return indexedBins;
}

void outputIndexedUniqueContinousBin(vecVecBin *indexedBins,vecBin *chromosomesContinousBin,int unitLen){
	size_t i,j,k;
	for(i=0;i<25;++i){
		for (j=0;j<indexedBins[i].n;++j) {
			printf("##%zd\t%zd\t%zd\n",j,unitLen*j+1,unitLen*(j+1));
			for (k=0;k<indexedBins[i].a[j].n;++k) {
				printf("%u\t%u\t%u\n",indexedBins[i].a[j].a[k].start,indexedBins[i].a[j].a[k].end,indexedBins[i].a[j].a[k].bin);
			}
			kv_destroy(indexedBins[i].a[j]);
		}
		kv_destroy(indexedBins[i]);
	}
	free(indexedBins);
}

void overlap(int *j,uint32_t chrLen ,int windows,vecBin *indexedBins,Bin *bin,uint16_t windowSize){
	uint32_t window_start,window_end;
	while (*j<=windows){
		window_start=windowSize*(*j)+1;
		window_end= (*j+1)*windowSize;
		window_end= (window_end>chrLen) ? chrLen : window_end;
		if (bin->end<window_start) {
			break;
		}
		else{
			if (bin->start<=window_start) {
				if (bin->end<=window_end) {
					if(indexedBins[*j].n){
						kv_push(Bin, indexedBins[*j], *bin);
					}else{
						kv_init(indexedBins[*j]);
						kv_push(Bin, indexedBins[*j], *bin);
					}
					break;
				}else{
					if(indexedBins[*j].n){
						kv_push(Bin, indexedBins[*j], *bin);
					}else{
						kv_init(indexedBins[*j]);
						kv_push(Bin, indexedBins[*j], *bin);
					}
					(*j)++;
				}
			}else{
				if (bin->start<=window_end) {
					if (bin->end<=window_end) {
						break;
					}else{
						if(indexedBins[*j].n){
							kv_push(Bin, indexedBins[*j], *bin);
						}else{
							kv_init(indexedBins[*j]);
							kv_push(Bin, indexedBins[*j], *bin);
						}
						(*j)++;
					}
				}else{
					(*j)++;
				}
			}
		}
	}
}

void IndexBins(dict *ht,uint16_t unitLen){
	if (ht){
		dictEntry** D=dump_dict(ht);
		unsigned long chri=0;
		int i;
		for (chri=0;chri<ht->used;chri++){
			ArrayBin *array=(ArrayBin *)(D[chri]->val);
			int window_count = array->bins[array->count-1].end/unitLen+1;
			vecBin *indexedBins = (vecBin *)calloc(window_count,sizeof(vecBin));
			int j=0;
			for (i=0;i<array->count ;i++ ) {
				overlap(&j,array->bins[array->count-1].end,window_count,indexedBins,array->bins+i,unitLen);
			}
			for (i=0;i<window_count;++i) {
				printf("##%d\t%d\t%d\n",i,unitLen*i+1,unitLen*(i+1));
				for (j=0;j<(int)indexedBins[i].n;++j){
					printf("%u\t%u\t%u\n",indexedBins[i].a[j].start,indexedBins[i].a[j].end,indexedBins[i].a[j].bin);
				}
				kv_destroy(indexedBins[i]);
			}
			free(indexedBins);
		}
		free(D);
	}
}
