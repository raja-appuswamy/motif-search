CPP=g++
CPPFLAGS=-g -Wall -pthread -std=c++17 
LIBS=-lz -lboost_program_options
ifneq ($(DEBUG),)
	CPPFLAGS+=-O0 -DDBGPRINT
else
	CPPFLAGS+=-O2
endif

HEADERS = $(wildcard *.h)
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
.PHONY = clean

all: motif-search

motif-search: $(OBJS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -o $@ $(OBJS) ${LIBS}

sample: motif-search
	./motif-search -r `ls sample_data/reads/EINS3-guppy2/*.fastq | xargs echo -n` -m sample_data/motifs/EINS.fa -l 5

clean:
	rm motif-search *.o
