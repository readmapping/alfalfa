CC=gcc
CXX=g++
FLAGS = -O3 -DNDEBUG -pthread
LIBS=-lm -lz
CCSRC = qsufsort.c ksw.c kthread.c
CXXSRC = main.cpp options.cpp sparseSA.cpp fasta.cpp dp.cpp utils.cpp mapper.cpp performanceUtils.cpp
PREFIX=/usr/local
src_dir = src
build_dir = build
sources += $(patsubst %.cpp, $(src_dir)/%.cpp, $(CXXSRC))
sources += $(patsubst %.c, $(src_dir)/%.c, $(CCSRC))
objs += $(patsubst %.cpp, $(build_dir)/%.cpp.o, $(CXXSRC))
objs += $(patsubst %.c, $(build_dir)/%.c.o, $(CCSRC))
EXEC = alfalfa

all: alfalfa

alfalfa: dir $(objs)
	$(CXX) $(FLAGS) $(objs) -o $@ $(LIBS)
	
dir:
	mkdir -p $(build_dir)

$(build_dir)/%.cpp.o: $(src_dir)/%.cpp
	$(CXX) $(FLAGS) -Wall -o $@ -c $<

$(build_dir)/%.c.o: $(src_dir)/%.c
	$(CC) $(FLAGS) -Wall -o $@ -c $<

install: alfalfa
	mkdir -p $(PREFIX)/bin
	cp alfalfa $(PREFIX)/bin/

# .PHONY assures clean is exected even if there is a file "./clean" in
# the directory. The same for doc.
.PHONY: clean doc
doc: 
	doxygen
clean: 
	rm -rf $(build_dir) 
	rm -f *~ .depend alfalfa

# Create all the dependencies between the source files. 
.depend:
	$(CXX) -MM $(sources) > .depend

# The - prevents make from complaining about a missing .depend
-include .depend
