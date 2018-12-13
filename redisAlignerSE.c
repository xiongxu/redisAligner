//gcc -g -O3 -Wall  klib/kthread.c redis_query_kthread_block_uniq.c -o redis_query_kthread_block_uniq -I./klib -I./hiredis -L./hiredis -lhiredis -lpthread -lz
#include <inttypes.h>
#include "hiredis.h"
#include "sds.h"
#include "uniqueBin2Hash.h"

#define N_PIPE 10000
#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define int2base(c) (1413956417 >> (8*(c))&0xff)
#define int2complbase(c) (1094928212 >> (8*(c))&0xff)

static uint8_t nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static uint8_t groups[256] = {
  1,7,7,0,2,10,6,8,3,11,8,5,10,8,6,6,
  5,0,4,1,2,7,8,11,4,8,3,8,1,7,11,3,
  4,3,4,4,7,6,6,7,1,5,1,6,0,5,9,8,
  0,7,5,5,3,2,10,9,10,9,10,11,11,3,7,9,
  7,3,5,11,3,2,6,8,8,11,10,7,3,8,6,7,
  5,4,6,3,6,0,1,0,2,10,0,10,6,9,5,3,
  1,5,2,6,10,6,9,5,3,2,2,1,10,9,6,11,
  0,9,2,1,11,0,2,6,4,1,10,10,6,6,11,7,
  9,11,0,9,8,11,6,6,7,5,2,0,6,6,8,3,
  6,10,10,4,4,9,9,5,0,9,2,1,11,2,1,11,
  1,7,3,4,8,5,7,1,11,11,9,10,8,9,2,5,
  5,2,2,2,2,10,4,3,9,5,6,3,9,10,9,2,
  2,11,6,7,5,9,2,1,1,5,1,8,10,9,10,4,
  3,10,7,10,2,10,11,2,9,8,2,5,4,11,10,8,
  1,4,9,2,11,0,8,3,2,0,8,0,3,5,3,4,
  0,5,8,1,5,0,8,9,4,6,2,8,11,8,0,4
};

static uint32_t centromere[24][3]={
  {121500000,128900000,249250621},
  {90500000,96800000,243199373},
  {87900000,93900000,198022430},
  {48200000,52700000,191154276},
  {46100000,50700000,180915260},
  {58700000,63300000,171115067},
  {58000000,61700000,159138663},
  {43100000,48100000,146364022},
  {47300000,50700000,141213431},
  {38000000,42300000,135534747},
  {51600000,55700000,135006516},
  {33300000,38200000,133851895},
  {16300000,19500000,115169878},
  {16100000,19100000,107349540},
  {15800000,20700000,102531392},
  {34600000,38600000,90354753},
  {22200000,25800000,81195210},
  {15400000,19000000,78077248},
  {24400000,28600000,59128983},
  {25600000,29400000,63025520},
  {10900000,14300000,48129895},
  {12200000,17900000,51304566},
  {58100000,63000000,155270560},
  {11600000,13400000,59373566},
};

typedef struct _SequenceAliagnmentResult_ {
  uint8_t chr ;
  uint8_t GC : 7;
  uint8_t strand : 1;
  uint32_t pos;
} SRA;

typedef kvec_t(SRA) vecSRA;
/*
typedef struct _RedisQueryKey_ {
  sds code;
  uint8_t DBi : 6;
  uint8_t read : 1;
  uint8_t strand : 1;
  SRA *sra;
  RedisQueryKey *mate;
} RedisQueryKey;
*/
typedef struct _RedisQueryKey_ {
  sds code;
  uint8_t DBi : 7;
  uint8_t strand : 1;
  SRA *sra;
} RedisQueryKey;

typedef struct {
  gzFile fp;
  uint32_t max_lines;
  int buf_size, n_threads;
  char *buf;
  redisContext **c;
  FILE *fout;
  uint64_t **RC,TRC;
  uint32_t **BRC;
  double **BGC;
  vecVecBin *indexedContinousBins;
  vecSRA sraArray;
} pipeline_t;

typedef struct {
  uint32_t n_lines;
  sds *lines;
  redisContext **c;
  RedisQueryKey *codeLines;
  uint32_t *groupCount;
  uint32_t *AccumulatedCount;
  vecSRA sraBlock;
} step_t;

struct globalArgs_t {
  char *infile1;
  char *outfile;
  int readCount;
  int stepCount,blockCount;
  char *continousBinFile;
} globalArgs;

#define revcompEncode(readLen,tmpsds,revcomp_packsds) \
do{ \
  int ii,endPos=readLen-1; \
  for(ii=endPos;ii>=0;--ii){ \
    _set_pac(revcomp_packsds,endPos-ii, 3-_get_pac(tmpsds, ii)); \
  } \
}while(0)

#define revcompEncode2(readLen,sequence,packsds) \
do{ \
  int ii,endPos=readLen-1; \
  packsds = sdsnewlen(NULL,readLen/4); \
  for(ii=endPos;ii>=0;--ii){ \
    _set_pac(packsds,endPos-ii, 3- nst_nt4_table[(uint8_t)sequence[ii]]); \
  } \
}while(0)

#define PrintSeqFrom_sds(seq,packSeq,readLen,f) \
do{ \
  int j=0; \
  for(;j<readLen;++j){ \
    seq[j]=int2base(_get_pac(packSeq,j)); \
  } \
  fprintf(f,"%s\n",seq); \
}while(0)

static inline int hasN(const char *a);
static inline int readNextNode(gzFile fq,char *buf,int buf_size,sds *line);
static inline int compCode(const void *a, const void *b);
static inline int compPos(const void *a, const void *b);
static inline void *worker_pipeline(void *shared, int step, void *in);
static inline void kthreadQuery(void *data,long i,int tid);
static inline void encodeDoubleStrand(void *_data, long i, int tid);

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

void display_usage(char * argv[]){
  char *buffer=(char* )malloc(10240*sizeof(char));
  const char* usage=
"\nCopyright (c) InfiniGenomics 2016-2017\n" \
"Contact: XiongXu <xux@infinigenomics.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for Spliting large text based gzip fastq files into smaller fixed sized ones.\n" \
"Usage: %s [-1 Infile1] [-2 Infile2] [-o OUTFILE] [-r reads_count_per_file] [-c continousBinFile] [-h] \n" \
"Example1:\n  %s -1 /home/xux/data/RS14_L8_I001.R1.clean.fastq.gz -r 500000 -o out -c uniqBin.bin2pos.gz\n" \
"\n" \
"   [-1 Infile1] = fastq Infile1.                              [required]\n" \
"   [-2 Infile2] = fastq Infile2.                              [required]\n" \
"   [-o OUTPUT] = OUTPUT file.default is \'out\'.              [option]\n" \
"   [-r Threshold] = reads_count_per_file.                     [option]\n" \
"   [-c continousBinFile] = continousBinFile.                  [option]\n" \
"   [-b blockCount] = blockCount,default is 12.                [option]\n" \
"   [-s stepCount] = stepCount,default is 1.                   [option]\n" \
"   [-d dedup or not] = whether to remove identical reads.     [option]\n" \
"   [-v ] = verbose (output mapping details).                  [option]\n" \
"   [-h]        = This helpful help screen.                    [option]\n" \
"\n";
  sprintf(buffer,usage,argv[0],argv[0]);
  fprintf(stderr,"%s",buffer);
  free(buffer);
  exit(1);
}

int hasN(const char *a) {
    for (;*a;++a) {
      if (nst_nt4_table[(uint8_t)*a] >3) return 1;
    }
    return 0;
}
/*
int readNextNode(gzFile fq,char *buf,int buf_size,sds *line) {
  int flag=0;
  do {
    buf=gzgets(fq,buf,buf_size);
    if (!gzeof(fq)) {
      buf=gzgets(fq,buf,buf_size);
      if(flag) sdsfree(*line);
      *line=sdsnewlen(buf,36);
      buf=gzgets(fq,buf,buf_size);
      buf=gzgets(fq,buf,buf_size);
      flag++;
    } else {
      return 0;
    }
  }while(hasN(*line) || strlen(buf) < 37);
  return 1;
}*/

int readNextNode(gzFile fq,char *buf,int buf_size,sds *line) {
  do {
    buf=gzgets(fq,buf,buf_size);
    if (!gzeof(fq)) {
      buf=gzgets(fq,buf,buf_size);
      *line=sdsnewlen(buf,36);
      buf=gzgets(fq,buf,buf_size);
      buf=gzgets(fq,buf,buf_size);
    } else {
      return 0;
    }
  }while(0);
  return 1;
}

int compCode(const void *a, const void *b) {
  return ((RedisQueryKey *)a)->DBi - ((RedisQueryKey *)b)->DBi;
}

int compPos(const void *a, const void *b) {
  SRA *c = (SRA *)a;
  SRA *d = (SRA *)b;
  return (((uint64_t)c->chr<<32|c->pos<<1|c->strand) > ((uint64_t)d->chr<<32|d->pos<<1|d->strand));
}

void encodeDoubleStrand(void *_data, long i, int tid) {// kt_for() callback
  step_t *step = (step_t*)_data;
  step->codeLines[i].code=sdsnewlen(NULL,9);
  int j;
  for (j=0;j<36;++j){
    _set_pac(step->codeLines[i].code,j,nst_nt4_table[(uint8_t)step->lines[i][j]]);
  }
  step->codeLines[i].DBi= groups[(uint8_t)step->codeLines[i].code[0]];
  step->codeLines[i].sra = step->sraBlock.a+i;
  step->codeLines[i].strand = 1;
  int revcompIndex = i + step->n_lines;
  step->codeLines[revcompIndex].code = sdsnewlen(NULL,9);
  revcompEncode(36,step->codeLines[i].code,step->codeLines[revcompIndex].code);
  if(!sdscmp(step->codeLines[revcompIndex].code,step->codeLines[i].code)){
    fprintf(stderr,"palindromic sequence: %s\n",step->lines[i]);
  }
  step->codeLines[revcompIndex].DBi= groups[(uint8_t)step->codeLines[revcompIndex].code[0]];
  step->codeLines[revcompIndex].sra = step->sraBlock.a+i;
  step->codeLines[revcompIndex].strand = 0;
  sdsfree(step->lines[i]);
}

void kthreadQuery(void *data,long i,int tid) {
  step_t *s=(step_t *)data;
  uint32_t j,k=0,m,curr,tmp=0;
  redisReply *reply;
  for(j=s->AccumulatedCount[i];k<s->groupCount[i];) {
    redisAppendCommand(s->c[i],"HGET %b %b",s->codeLines[j].code,3,s->codeLines[j].code+3,6);
    ++k;++j;
    if(k % N_PIPE ==0 || k==s->groupCount[i] ) {
      uint32_t rest= k % N_PIPE ==0 &&  k!=s->groupCount[i] ? N_PIPE:k % N_PIPE;
      for(m=0;m<rest;++m){
        redisGetReply(s->c[i],(void **)&reply);
        curr = j-rest+m;
        if (reply && reply->type == REDIS_REPLY_STRING) {
          uint8_t str0=(uint8_t)reply->str[0],
                  str1=(uint8_t)reply->str[1];
          s->codeLines[curr].sra->chr = str0 & 0x1f ;
          s->codeLines[curr].sra->GC = (str0>>1&0x70) | (str1>>4);
          s->codeLines[curr].sra->strand = s->codeLines[curr].strand;
          s->codeLines[curr].sra->pos = (uint32_t)(str1 & 0xf);
          for (tmp=2;tmp<5 ;tmp++){
            s->codeLines[curr].sra->pos = s->codeLines[curr].sra->pos << 8 | (uint8_t)reply->str[tmp];
          }
        }
        sdsfree(s->codeLines[curr].code);
        freeReplyObject(reply);
      }
    }
  }
}

void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
  pipeline_t *p = (pipeline_t*)shared;
  uint32_t i=0;
  if (step == 0) { // step 0: read lines into the buffer
    step_t *s;
    s = calloc(1, sizeof(step_t));
    s->n_lines=0;
    s->lines = calloc(p->max_lines, sizeof(sds));
    s->c = p->c;
    s->groupCount = calloc(p->n_threads, sizeof(int));
    s->AccumulatedCount = calloc(p->n_threads, sizeof(int));
    while (readNextNode(p->fp,p->buf,p->buf_size,s->lines+s->n_lines)){
      if(hasN(*(s->lines+s->n_lines)) || strlen(p->buf) < 37){
        sdsfree(*(s->lines+s->n_lines));
        continue;
      }
      if (++s->n_lines >= p->max_lines) break;
    }
    if (s->n_lines) {
      s->codeLines = calloc(s->n_lines *2 , sizeof(RedisQueryKey));
      s->sraBlock.m=s->n_lines+1;
      s->sraBlock.n=0;
      s->sraBlock.a= (SRA *)calloc(s->sraBlock.m,sizeof(SRA));
      return s;
    }
  } else if (step == 1) {
    step_t *s = (step_t*)in;
    kt_for(p->n_threads, encodeDoubleStrand, in, s->n_lines);
    free(s->lines);
    return s;
  } else if (step == 2) {
    step_t *s = (step_t*)in;
    uint32_t k=s->n_lines*2,sumCount=0;
    for (i=0;i<k;++i) {
      s->groupCount[s->codeLines[i].DBi]++;
    }
    for (i=0;i<(uint32_t)p->n_threads;++i){
      s->AccumulatedCount[i]=sumCount;
      sumCount+=s->groupCount[i];
    }
    return s;
  } else if (step == 3) {
    qsort(((step_t*)in)->codeLines,((step_t*)in)->n_lines*2,sizeof(RedisQueryKey),compCode);
    return in;
  } else if (step == 4) {
    kt_for(p->n_threads, kthreadQuery ,in, p->n_threads);
    return in;
  } else if (step == 5) {
    step_t *s = (step_t*)in;
    if (s->n_lines) p->TRC += s->n_lines;
    for (i=0;i<s->n_lines;++i) {
      if (s->sraBlock.a[i].chr){
        kv_push(SRA, p->sraArray, s->sraBlock.a[i]);
      }
    }
    free(s->codeLines);
    free(s->groupCount);
    free(s->AccumulatedCount);
    kv_destroy(s->sraBlock);
    free(s);
  }
  return 0;
}

int main(int argc, char *argv[]) {
  long long begin=usec();
  int opt = 0, verbose=0, dedup=0;;
  globalArgs.infile1=NULL;
  globalArgs.outfile="out";
  globalArgs.readCount=500000;
  globalArgs.blockCount=12;
  globalArgs.stepCount=1;
  globalArgs.continousBinFile=NULL;
  const char *optString = "1:o:r:c:b:s:vdh?";
  if (argc<2) display_usage(argv);
  opt = getopt( argc, argv, optString);
  while( opt != -1 ) {
    switch( opt ) {
      case '1':
        globalArgs.infile1 = optarg;
        break;
      case 'c':
        globalArgs.continousBinFile = optarg;
        break;
      case 'o':
        globalArgs.outfile = optarg;
        break;
      case 'r':
        globalArgs.readCount = atoi(optarg);
        break;
      case 'b':
        globalArgs.blockCount = atoi(optarg);
        break;
      case 's':
        globalArgs.stepCount = atoi(optarg);
        break;
      case 'v':
        verbose++;
        break;
      case 'd':
        dedup++;
        break;
      case '?': /* fall-through is intentional */
      case 'h':
        display_usage(argv);
        break;
      default:
        fprintf(stderr,"error parameter!\n");
        break;
    }
    opt = getopt( argc, argv, optString );
  }
  pipeline_t pl;
  pl.fp=open_input_stream(globalArgs.infile1);
  int pl_threads =6,steps=6;
  uint32_t i=0,j=0,k=0;
  pl.buf_size = 0x200;
  pl.max_lines = globalArgs.readCount;
  pl.buf = calloc(pl.buf_size, 1);
  pl.n_threads = 12;
  pl.TRC=0;
  kv_init(pl.sraArray);
  pl.RC = calloc(24,sizeof(uint64_t *));
  pl.BRC = calloc(24,sizeof(uint32_t *));
  pl.BGC = calloc(24,sizeof(double *));
  vecBin *continousBin = loadUnique2BinFile(globalArgs.continousBinFile);
  pl.indexedContinousBins = indexUniqueContinousBin(continousBin,1024);
  uint32_t maxCount = 0;
  for (i=0;i<24;++i){
    if(maxCount < continousBin[i].n ) maxCount = continousBin[i].n;
  }
  for (i=0;i<24;++i) {
    pl.BRC[i] = (uint32_t *)calloc(maxCount,sizeof(uint32_t));
    pl.BGC[i] = (double *)calloc(maxCount,sizeof(double));
    pl.RC[i] = (uint64_t *)calloc(6,sizeof(uint64_t));
  }
  pl.c = (redisContext **)calloc(pl.n_threads,sizeof(redisContext *));
  for (i=0;i<(uint32_t)pl.n_threads;++i){
    pl.c[i] = redisConnect((char *)"127.0.0.1", 6379+i);
    if (pl.c[i]->err) error("Error: %s\n", pl.c[i]->errstr);
  }
  pl.fout=fcreat_outfile(globalArgs.outfile,".redis.tsv");
  kt_pipeline(pl_threads, worker_pipeline, &pl, steps);
  if (dedup){
    qsort(pl.sraArray.a,pl.sraArray.n,sizeof(SRA),compPos);
    SRA l_sra={0};
    uint64_t redundency=0;
    for (i=0;i<pl.sraArray.n;++i){
      if( ((uint64_t)l_sra.chr<<32|l_sra.pos<<1|l_sra.strand) == ((uint64_t)pl.sraArray.a[i].chr<<32|pl.sraArray.a[i].pos<<1|pl.sraArray.a[i].strand) ){
        ++redundency;
      }else{
        uint8_t chrIndex=pl.sraArray.a[i].chr-1;
        pl.RC[chrIndex][0]++;
        pl.RC[chrIndex][1]+=pl.sraArray.a[i].GC;
        if (pl.sraArray.a[i].pos < centromere[chrIndex][0]){
          pl.RC[chrIndex][2]++;
          pl.RC[chrIndex][3]+=pl.sraArray.a[i].GC;
        }else if(pl.sraArray.a[i].pos > centromere[chrIndex][1]) {
          pl.RC[chrIndex][4]++;
          pl.RC[chrIndex][5]+=pl.sraArray.a[i].GC;
        }
        uint32_t binIndex = pl.sraArray.a[i].pos/1024,ii=0;
        for (ii=0;ii<pl.indexedContinousBins[chrIndex].a[binIndex].n;++ii){
          if (pl.sraArray.a[i].pos > pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].start &&
            pl.sraArray.a[i].pos <= pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].end){
            uint32_t CurrentBlock = pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].bin;
            pl.BRC[chrIndex][CurrentBlock]++;
            pl.BGC[chrIndex][CurrentBlock]+=pl.sraArray.a[i].GC;
            break;
          }
        }
      }
      l_sra = pl.sraArray.a[i];
    }
    fprintf(stderr,"redundency ratio: %" PRIu64 " / %" PRIu64 " = %%%lf\n", redundency,pl.TRC,1.0*redundency/pl.TRC*100);
  }else{
    for (i=0;i<pl.sraArray.n;++i){
      uint8_t chrIndex=pl.sraArray.a[i].chr-1;
      pl.RC[chrIndex][0]++;
      pl.RC[chrIndex][1]+=pl.sraArray.a[i].GC;
      if (pl.sraArray.a[i].pos < centromere[chrIndex][0]){
        pl.RC[chrIndex][2]++;
        pl.RC[chrIndex][3]+=pl.sraArray.a[i].GC;
      }else if(pl.sraArray.a[i].pos > centromere[chrIndex][1]) {
        pl.RC[chrIndex][4]++;
        pl.RC[chrIndex][5]+=pl.sraArray.a[i].GC;
      }
      uint32_t binIndex = pl.sraArray.a[i].pos/1024,ii=0;
      for (ii=0;ii<pl.indexedContinousBins[chrIndex].a[binIndex].n;++ii){
        if (pl.sraArray.a[i].pos > pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].start &&
          pl.sraArray.a[i].pos <= pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].end){
          uint32_t CurrentBlock = pl.indexedContinousBins[chrIndex].a[binIndex].a[ii].bin;
          pl.BRC[chrIndex][CurrentBlock]++;
          pl.BGC[chrIndex][CurrentBlock]+=pl.sraArray.a[i].GC;
          break;
        }
      }
    }
  }
  if (verbose){
    for (i=0;i<pl.sraArray.n;++i){
      fprintf(stdout,"%" PRIu8 "\t%" PRIu32 "\t%" PRIu8 "\t%c\n",pl.sraArray.a[i].chr,pl.sraArray.a[i].pos,pl.sraArray.a[i].GC,pl.sraArray.a[i].strand?'+':'-');
    }
  }
  kv_destroy(pl.sraArray);

  uint64_t noXYRC=0,noXYRCpter=0,noXYRCqter=0;
  for (i=0;i<22;++i) {
    noXYRC+=pl.RC[i][0];
    noXYRCpter+=pl.RC[i][2];
    noXYRCqter+=pl.RC[i][4];
  }
  uint64_t sumRC=noXYRC;
  for (;i<24;++i){ sumRC+=pl.RC[i][0];}
  double chrYPCT = (double)pl.RC[23][0]/noXYRC;
  FILE *blkOut = fcreat_outfile(globalArgs.outfile,".blk");
  //fprintf(blkOut,"#%s\n",chrYPCT<0.00007?"F":"M");
  fprintf(pl.fout,"chr\tNRC\tRCchr\tGCchr\tNRCpter\tRCpter\tGCpter\tNRCqter\tRCqter\tGCqter\n");
  for (i=0;i<24;++i){
    fprintf(pl.fout,"%" PRIu32 "\t%lf\t%" PRIu64 "\t%lf\t%lf\t%" PRIu64 "\t%lf\t%lf\t%" PRIu64 "\t%f\n",i+1,
      (double)pl.RC[i][0]/noXYRC,pl.RC[i][0],(double)pl.RC[i][1]/(36*pl.RC[i][0]),
      (double)pl.RC[i][2]/noXYRCpter,pl.RC[i][2],pl.RC[i][2] ? (double)pl.RC[i][3]/(36*pl.RC[i][2]) : 0,
      (double)pl.RC[i][4]/noXYRCqter,pl.RC[i][4],pl.RC[i][4] ? (double)pl.RC[i][5]/(36*pl.RC[i][4]) : 0);
    fprintf(blkOut, "%" PRIu32 "\t%zu",i+1,continousBin[i].n);
    for (j=0;j<maxCount;j+=globalArgs.stepCount){
      //fprintf(blkOut, "\t%" PRIu32 "\t%lf",pl.BRC[i][j],pl.BRC[i][j]==0? 0: pl.BGC[i][j]/(pl.BRC[i][j]*36.0)*100);
      uint32_t BsumRC=0;
      double BsumGC=0;
      uint32_t maxK = j+globalArgs.blockCount <=maxCount ? j+globalArgs.blockCount : maxCount;
      for (k=j;k<maxK;++k){
        BsumRC+= pl.BRC[i][k];
        BsumGC+= pl.BGC[i][k];
      }
      fprintf(blkOut, "\t%lf\t%lf",1.0*BsumRC/noXYRC*1000000,BsumRC==0 ? 0: BsumGC/(BsumRC*36.0)*100);
    }
    fputc('\n',blkOut);
  }
  fclose(blkOut);
  fprintf(stderr,"mapping ratio: %" PRIu64 " / %" PRIu64 " = %%%lf\n", sumRC,pl.TRC,1.0*sumRC/pl.TRC*100);
  fprintf(stderr,"chrY%%: %" PRIu64 " / %" PRIu64 " = %lf, SexType: %s\n",pl.RC[23][0],noXYRC,chrYPCT,chrYPCT<0.00007?"F":"M");
  for (i=0;i<24;++i) {
    free(pl.RC[i]);
    free(pl.BRC[i]);
    free(pl.BGC[i]);
  }
  free(pl.RC);
  free(pl.BRC);
  free(pl.BGC);
  for (i=0;i<(uint32_t)pl.n_threads;++i) redisFree(pl.c[i]);
  free(pl.buf);
  gzclose(pl.fp);
  if (pl.fout!=stdout) fclose(pl.fout);
  fprintf(stderr,"finished in : %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
  return 0;
}
