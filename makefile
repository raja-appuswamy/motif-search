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
.PHONY = depend clean

all: motif-search

motif-search: $(OBJS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -o $@ $(OBJS) ${LIBS}

clean:
	rm motif-search *.o
