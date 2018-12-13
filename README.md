# redisAligner
NGS 36-mer perfect match sequecing aligner based on redis

## Installation
1. requirement
```
yum -y install libbz2-devel liblzma-devel ncurses-devel openssl-devel libcurl-devel 
wget ftp://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph;
chmod 750 bigWigToBedGraph;
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz;
tar -zxvf bedtools-2.27.1.tar.gz;
cd bedtools2;
make;
perl -e 'map { print "wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr$_.fa.gz && gzip -cd chr$_.fa.gz > chr$_.fa\n"} (1..22,"X","Y","M")'|bash;
perl -e 'print "cat ",join(" ",map {"chr$_.fa" }(1..22,"X","Y","M"))," > hg19.fa\n";'|bash;
```
2. compilation
```
make
```

## build index of hg19
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig;
bigWigToBedGraph wgEncodeCrgMapabilityAlign36mer.bigWig wgEncodeCrgMapabilityAlign36mer.bedGraph;
bedtools sort -faidx <(perl -e 'print join("\n",map {"chr".$_} (1..22,"X","Y","M")),"\n";') -i <(awk -F"\t" '($4==1){print $0}' wgEncodeCrgMapabilityAlign36mer.bedGraph) > wgEncodeCrgMapabilityAlign36mer_1.sort.BedGraph
#this step may cost 221G memory in 10 hours
./ref2aof -i wgEncodeCrgMapabilityAlign36mer_1.BedGraph -o hg19 -r hg19.fasta 
```
## load data to redis
1. copy all the generated aof files into redis-4.0.9 directory
2. excute sh batch_start_aof.sh, if you have generated rdb files,you just need to excute sh batch_start_rdb.sh, it is faster than load aof file into redis.

## query fastq file
```
./redis_query_kthread_block_uniqPE -1 file_R1.fastq.gz -2 file_R2.fastq.gz -o file_out -c uniqBin.bin2pos.gz -d  2> file.log
```
