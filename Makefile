CC=			gcc
CFLAGS=		-O3 -Wall -W -Wstrict-prototypes -Wwrite-strings -g -ggdb
LDFLAGS=	-lhiredis -lz -lpthread
INCLUDES=	-I$(CURDIR)/hiredis -I$(CURDIR)/zlib-1.2.11/include -I$(CURDIR)/klib -I.
LIBPATH=    -L$(CURDIR)/hiredis -L$(CURDIR)/zlib-1.2.11/lib -L.

.PHONY:all clean
all: redisAlignerSE redisAlignerPE ref2aof uniqueBin2Hash
#SUBDIRS = `find -maxdepth 1 -type d | sed "1d"`
$(CURDIR)/hiredis/libhiredis.a:
	wdir=`pwd`; \
	cd  $(CURDIR)/hiredis;\
	$(MAKE)  ;\
	cd $$wdir

$(CURDIR)/zlib-1.2.11:
	wget http://www.zlib.net/zlib-1.2.11.tar.gz; \
	tar -zxvf zlib-1.2.11.tar.gz && mv zlib-1.2.11 zlib;\
	cd zlib;\
	test -d $(CURDIR)/zlib-1.2.11 || mkdir -p $(CURDIR)/zlib-1.2.11;\
	./configure --prefix=$(CURDIR)/zlib-1.2.11;\
	make;\
	make install;\
	cd $(CURDIR) && rm -rf zlib

$(CURDIR)/htslib/libhts.a: $(CURDIR)/zlib-1.2.11 $(CURDIR)/htslib
	wdir=`pwd`; \
	cd  $(CURDIR)/htslib;\
	$(MAKE) INCLUDES="-I. -I$(CURDIR)/zlib-1.2.11/include";\
	cd $$wdir

redisAlignerSE:redisAlignerSE.c klib/kthread.c $(CURDIR)/zlib-1.2.11 $(CURDIR)/hiredis/libhiredis.a
	$(CC) $(CFLAGS) klib/kthread.c $< -o $@ $(INCLUDES) $(LIBPATH) $(LDFLAGS)

redisAlignerPE:redisAlignerPE.c klib/kthread.c $(CURDIR)/zlib-1.2.11 $(CURDIR)/hiredis/libhiredis.a
	$(CC) $(CFLAGS) klib/kthread.c $< -o $@ $(INCLUDES) $(LIBPATH) $(LDFLAGS) -lm

ref2aof:ref2aof.c $(CURDIR)/zlib-1.2.11 $(CURDIR)/hiredis/libhiredis.a $(CURDIR)/htslib/libhts.a
	$(CC) $(CFLAGS) $< -o $@ $(INCLUDES) -I$(CURDIR)/htslib/htslib $(LIBPATH) -L$(CURDIR)/htslib $(LDFLAGS) -lhts

uniqueBin2Hash:uniqueBin2Hash.c $(CURDIR)/hiredis/libhiredis.a $(CURDIR)/zlib-1.2.11
	$(CC) $(CFLAGS) $< -o $@ $(INCLUDES) $(LIBPATH) $(LDFLAGS)

clean:
	rm -rf zlib-1.2.11 zlib-1.2.11.tar.gz *.o redisAlignerSE redisAlignerPE ref2aof uniqueBin2Hash *.dSYM;\
	wdir=`pwd`; \
	cd ./hiredis;\
	$(MAKE) clean;\
	cd $$wdir;\
	cd  $(CURDIR)/htslib;\
	$(MAKE) clean ;\
	cd $$wdir
