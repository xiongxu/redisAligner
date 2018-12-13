//gcc -g -static -O3 -Wall uniqueBin2Hash.c -o uniqueBin2Hash -I./klib -I./hiredis -L./hiredis -lhiredis -lz
#include <getopt.h>
#include <inttypes.h>
#include "uniqueBin2Hash.h"

struct globalArgs_t {
  char *infile1;
  char *outfile;
  int windowSize;
} globalArgs;

void display_usage(char * argv[]){
  char *buffer=(char* )malloc(10240*sizeof(char));
  const char* usage=
"\nCopyright (c) 2017\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for generating uniqBin.bin2pos File.\n" \
"Usage: %s [-1 unique mappability File] [-w Bin Size] [-o OUTFILE] [-h] \n" \
" gzip -cd wgEncodeCrgMapabilityAlign36mer_1.BedGraph.gz | awk -F\"\\t\" 'BEGIN{OFS=\"\\t\"}{if(chr!=$1){sumLen=0}chr=$1;Len=$3-$2;sumLen+=Len;printf(\"%%s\\t%%d\\t%%d\\t%%d\\t%%d\\t%%d\\t%%d\\n\",$1,$2,$3,$4,Len,sumLen,sumLen/5000)}' | gzip -c > wgEncodeCrgMapabilityAlign36mer_2.BedGraph.gz \n" \
"Example1:\n  %s -1 wgEncodeCrgMapabilityAlign36mer_1.BedGraph.gz -o uniqBin.bin2pos  > uniqBin.bin2pos.index -w 5000 \n" \
"\n" \
"   [-1 READ1]  = unique mappability File.                             [required]\n" \
"   [-w READ2]  = Bin Size.                                            [option]\n" \
"   [-o OUTPUT] = OUTPUT file.                                         [required]\n" \
"   [-h]        = This helpful help screen.                            [option]\n" \
"\n";
  sprintf(buffer,usage,argv[0],argv[0]);
  fprintf(stderr,"%s",buffer);
  free(buffer);
  exit(1);
}

int main(int argc, char *argv[]) {
  int opt = 0;
  globalArgs.infile1=NULL;
  globalArgs.outfile="out";
  globalArgs.windowSize=5000;
  const char *optString = "1:o:w:h?";
  if (argc<2) display_usage(argv);
  opt = getopt( argc, argv, optString);
  while( opt != -1 ) {
    switch( opt ) {
      case '1':
        globalArgs.infile1 = optarg;
        break;
      case 'o':
        globalArgs.outfile = optarg;
        break;
      case 'w':
        globalArgs.windowSize = atoi(optarg);
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
  dictType *mydictType=initDictType();
  dict *ht = load_bin(globalArgs.infile1,mydictType,globalArgs.windowSize);
  outputSortedHash(ht,globalArgs.outfile);
  IndexBins(ht,1024);
  // vecBin *continousBin = loadUnique2BinFile(globalArgs.infile1);
  // vecVecBin *indexedContinousBins = indexUniqueContinousBin(continousBin,1024);
  // outputIndexedUniqueContinousBin(indexedContinousBins,continousBin,1024);
  return 0;
}
