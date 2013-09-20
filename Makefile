#MAC=1
all		: cortex_con
#all	: cortex_con_cs solid_converter 
#all: solid_converter
ifndef CC
  CC = gcc	
endif

ifdef MAC
  CC = gcc
endif

BIN = bin
LIB = lib

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

ifeq ($(MAXK),160)
   BITFIELDS = 5
endif

ifeq ($(MAXK),192)
   BITFIELDS = 6
endif

ifeq ($(MAXK),223)
   BITFIELDS = 7
endif

ifeq ($(MAXK),255)
   BITFIELDS = 8
endif

ifndef BITFIELDS
   BITFIELDS = 1
   MAXK = 31
endif 

# Main program includes
IDIR_BASIC =include/basic
IDIR_UTIL =include/util
IDIR_HASH  =include/hash_table 
IDIR_CORTEX_CON = include/cortex_con
IDIR_ALIGNMENT = include/alignment
IDIR_STATS = include/stats

# Test code includes
IDIR_BASIC_TESTS =include/test/basic
IDIR_HASH_TABLE_TESTS =include/test/hash_table
IDIR_CORTEX_CON_TESTS=include/test/cortex_con

#Default CUnit installation path in some Linux distributions. 
IDIR_CUNIT=/usr/local/include/CUnit/ 


UNAME := $(shell uname)


ifeq ($(UNAME), Darwin)
    MAC=1
endif

# Mac OS X specific. Assuming CUnit is installed with MacPorts
ifdef MAC
MACFLAG = -fnested-functions 
#-L/opt/local/lib/ 
IDIR_CUNIT=/opt/local/include/CUnit/ 
CFLAGS_CUNIT = -L/opt/local/lib/ -lncurses
#LD_CUNIT=
endif 

# 64bit architecture by default
ARCH =  -m64 
ifdef 32_BITS
 ARCH =  
endif

# Compiler options
OPT		= $(ARCH) -Wall -O3 $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g

ifdef DEBUG
OPT	= $(ARCH) -Wall -O0 $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -g -pthread
endif

ifdef ENABLE_READ_PAIR
 OPT +=  -DENABLE_READ_PAIR
endif

ifdef ENABLE_READ_PAIR_OLD
 OPT += -DENABLE_READ_PAIR_OLD
endif

ifdef READ_PAIR_DEBUG_GRAPH
 OPT += -DREAD_PAIR_DEBUG_GRAPH
endif

ifdef DEBUG_PRINT_LABELS
 OPT += -DDEBUG_PRINT_LABELS
endif

metacortex:             OPT += -DNUMBER_OF_COLOURS=1 -DMETACORTEX
cortex_con:		OPT += -DNUMBER_OF_COLOURS=1 
cortex_con_rp:		OPT += -DNUMBER_OF_COLOURS=1 -DENABLE_READ_PAIR
cortex_con_mp:		OPT += -DNUMBER_OF_COLOURS=1 -DENABLE_MARK_PAIR
cortex_bub:		OPT += -DNUMBER_OF_COLOURS=2 -DENABLE_BUBBLEPARSE
bubbleparse:		OPT += -DNUMBER_OF_COLOURS=2 -DENABLE_BUBBLEPARSE -DINCLUDE_QUALITY_SCORES
graphout:		OPT += -DNUMBER_OF_COLOURS=2
kmerinfo:		OPT += -DNUMBER_OF_COLOURS=2
solid_converter: OPT += -DSOLID
cortex_con_cs:	OPT += -DSOLID
kmer_contamination:	OPT += -DKMER_CONTAMINATION -DNUMBER_OF_COLOURS=2 -DKMER_TOOLS
kmer_hash_build:	OPT += -DKMER_HASH_BUILD -DKMER_TOOLS	
kmer_filter:	OPT += -DKMER_TOOLS	
# Include dirs
CFLAGS_CORTEX_CON	= -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_UTIL)
CFLAGS_CORTEX_BUB	= -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_UTIL)
CFLAGS_GRAPHOUT         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_FILTERREADS         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_BUBBLEPARSE      = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH)
CFLAGS_KMERINFO         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH)
CFLAGS_LOGPARSE         = -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_HASH_TABLE_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH_TABLE_TESTS) 
CFLAGS_CORTEX_CON_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_CORTEX_CON_TESTS) 
CFLAGS_BASIC_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_BASIC_TESTS)
CFLAGS_KMER_STATS   	= -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL) -I$(IDIR_STATS)
CFLAGS_KMER_HASH_BUILD   	= -I$(IDIR_BASIC) -I$(IDIR_CORTEX_CON) -I$(IDIR_HASH) -I$(IDIR_UTIL) -I$(IDIR_STATS)
CFLAGS_DEMULTIPLEXER	= -I$(IDIR_BASIC) -I$(IDIR_UTIL) -I$(IDIR_STATS)
# Program objects
CORTEX_CON_OBJ = obj/cortex_con/file_format.o obj/cortex_con/flags.o obj/cortex_con/cleaning.o obj/cortex_con/path.o obj/cortex_con/perfect_path.o obj/cortex_con/branches.o obj/cortex_con/y_walk.o obj/cortex_con/cmd_line.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/cortex_con/cortex_con.o obj/cortex_con/logger.o obj/cortex_con/metacortex.o obj/cortex_con/coverage_walk.o obj/util/node_queue.o

ifdef ENABLE_READ_PAIR
CORTEX_CON_OBJ += obj/cortex_con/binary_tree.o obj/cortex_con/read_pair.o
BIN_SUFFIX = $(join rp_,$(MAXK))
else
BIN_SUFFIX = $(MAXK)
endif 
cortex_con_rp:		CORTEX_CON_OBJ +=  obj/cortex_con/binary_tree.o obj/cortex_con/read_pair.o
cortex_con_mp:		CORTEX_CON_OBJ +=  obj/cortex_con/mark_pair.o
CORTEX_BUB_OBJ = $(CORTEX_CON_OBJ) obj/cortex_con/bubble_find.o

BUBBLEPARSE_OBJ = obj/cortex_con/flags.o obj/cortex_con/path.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_table.o obj/cortex_con/file_reader.o obj/util/bubbleparse.o obj/cortex_con/logger.o obj/cortex_con/hash_value.o obj/cortex_con/file_format.o
KMERINFO_OBJ =  obj/cortex_con/flags.o obj/cortex_con/path.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_table.o obj/cortex_con/file_reader.o obj/util/kmerinfo.o obj/cortex_con/logger.o obj/cortex_con/hash_value.o
GRAPHOUT_OBJ = obj/cortex_con/flags.o obj/cortex_con/path.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_table.o obj/cortex_con/file_reader.o obj/cortex_con/dB_graph.o obj/util/graphout.o obj/cortex_con/perfect_path.o obj/cortex_con/logger.o obj/cortex_con/hash_value.o obj/util/graph_formats.o obj/util/node_queue.o obj/cortex_con/cleaning.o obj/cortex_con/coverage_walk.o obj/util/graph_tools.o obj/cortex_con/file_format.o
FILTERREADS_OBJ = obj/util/filter_reads.o obj/cortex_con/flags.o obj/cortex_con/path.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_table.o obj/cortex_con/file_reader.o obj/cortex_con/dB_graph.o obj/cortex_con/perfect_path.o obj/cortex_con/logger.o obj/cortex_con/hash_value.o obj/util/node_queue.o obj/cortex_con/cleaning.o obj/cortex_con/file_format.o
LOGPARSE_OBJ = obj/util/cortex_log_parse.o obj/util/reformat_log.o
BASIC_TESTS_OBJ	= obj/test/binary_kmer.o obj/test/seq.o obj/test/test_binary_kmer.o obj/test/test_seq.o obj/test/run_basic_tests.o  obj/cortex_con/logger.o
HASH_TABLE_TESTS_OBJ = obj/cortex_con/flags.o obj/test/run_hash_table_tests.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/test/test_hash.o obj/cortex_con/binary_kmer.o  obj/cortex_con/seq.o obj/cortex_con/logger.o
GRAPH_TESTS_OBJ = obj/cortex_con/branches.o obj/cortex_con/file_format.o obj/test/test_dB_graph.o obj/cortex_con/logger.o  obj/cortex_con/cleaning.o  obj/cortex_con/perfect_path.o obj/cortex_con/path.o obj/cortex_con/flags.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/cortex_con/y_walk.o  obj/test/test_file_reader.o obj/test/test_graph_element.o obj/test/run_dB_graph_tests.o 
SOLID_CONV_OBJ = obj/cortex_con/seq.o obj/cortex_con/colour_space_seq.o obj/cortex_con/binary_kmer.o obj/cortex_con/logger.o obj/util/solid_converter.o
KMER_COUNT_OBJ = obj/cortex_con/flags.o obj/cortex_con/cleaning.o obj/cortex_con/path.o obj/cortex_con/perfect_path.o  obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/stats/kmer_stats.o obj/stats/kmer_hash.o obj/cortex_con/logger.o obj/stats/kmer_reader.o obj/cortex_con/file_format.o
KMER_HASH_BUILD_OBJ = obj/cortex_con/flags.o obj/cortex_con/file_format.o obj/cortex_con/cleaning.o obj/cortex_con/path.o obj/cortex_con/perfect_path.o  obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/stats/kmer_hash.o obj/stats/kmer_hash_build.o obj/stats/kmer_reader.o  obj/cortex_con/logger.o
KMER_FILTER_OBJ = obj/cortex_con/seq_io.o obj/cortex_con/flags.o obj/cortex_con/file_format.o obj/cortex_con/cleaning.o obj/cortex_con/path.o obj/cortex_con/perfect_path.o  obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/stats/kmer_hash.o obj/stats/kmer_filter.o obj/stats/kmer_reader.o  obj/cortex_con/logger.o
SUBSAMPLER_OBJ = obj/cortex_con/flags.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/logger.o obj/cortex_con/file_format.o obj/util/subsampler.o
DEMULTIPLEXER_OBJ = obj/cortex_con/file_format.o  obj/cortex_con/flags.o obj/cortex_con/seq_io.o obj/cortex_con/seq.o obj/cortex_con/logger.o obj/cortex_con/file_format.o obj/stats/demultiplexer.o

#Library objects
LIBRARY_OBJ =  obj/cortex_con/file_format.o obj/cortex_con/analysis.o obj/cortex_con/flags.o obj/cortex_con/cleaning.o obj/cortex_con/path.o obj/cortex_con/perfect_path.o obj/cortex_con/branches.o obj/cortex_con/y_walk.o obj/cortex_con/cmd_line.o obj/cortex_con/binary_kmer.o obj/cortex_con/seq.o obj/cortex_con/element.o obj/cortex_con/hash_value.o obj/cortex_con/hash_table.o obj/cortex_con/dB_graph.o obj/cortex_con/file_reader.o obj/cortex_con/cortex_con.o obj/cortex_con/logger.o obj/cortex_con/metacortex.o obj/cortex_con/coverage_walk.o obj/util/node_queue.o

# Main rules
cortex_con_cs : remove_objects $(CORTEX_CON_OBJ) 
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_con_cs_$(BIN_SUFFIX) $(CORTEX_CON_OBJ)

solid_converter: remove_objects $(SOLID_CONV_OBJ) 
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/solid_converter $(SOLID_CONV_OBJ)

metacortex : remove_objects $(CORTEX_CON_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/metacortex_k$(BIN_SUFFIX) $(CORTEX_CON_OBJ)

cortex_con : remove_objects $(CORTEX_CON_OBJ) 
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_con_$(BIN_SUFFIX) $(CORTEX_CON_OBJ)

cortex_con_rp : remove_objects $(CORTEX_CON_OBJ)  obj/cortex_con/binary_tree.o obj/cortex_con/read_pair.o
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_con_rp_$(BIN_SUFFIX) $(CORTEX_CON_OBJ)

cortex_con_mp : remove_objects $(CORTEX_CON_OBJ)  obj/cortex_con/mark_pair.o
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_con_mp_$(BIN_SUFFIX) $(CORTEX_CON_OBJ)

cortex_bub : remove_objects $(CORTEX_BUB_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_bub_$(BIN_SUFFIX) $(CORTEX_BUB_OBJ)

kmerinfo: remove_objects $(KMERINFO_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmerinfo $(KMERINFO_OBJ)

graphout:remove_objects $(GRAPHOUT_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/graphout_$(BIN_SUFFIX) $(GRAPHOUT_OBJ)

filterreads:remove_objects $(FILTERREADS_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/filterreads_$(BIN_SUFFIX) $(FILTERREADS_OBJ)

bubbleparse: remove_objects $(BUBBLEPARSE_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/bubbleparse_$(MAXK) $(BUBBLEPARSE_OBJ)

subsampler: remove_objects $(SUBSAMPLER_OBJ)
		mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/subsampler $(SUBSAMPLER_OBJ)

demultiplexer: remove_objects $(DEMULTIPLEXER_OBJ)
		mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/demultiplexer $(DEMULTIPLEXER_OBJ)

reformat_log: remove_objects $(LOGPARSE_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/reformat_log $(LOGPARSE_OBJ)

kmer_contamination: remove_objects $(KMER_COUNT_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_contamination_$(BIN_SUFFIX) $(KMER_COUNT_OBJ)

kmer_hash_build: remove_objects $(KMER_HASH_BUILD_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_hash_build_$(BIN_SUFFIX) $(KMER_HASH_BUILD_OBJ)

kmer_filter: remove_objects $(KMER_FILTER_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/kmer_filter_$(BIN_SUFFIX) $(KMER_FILTER_OBJ)

run_basic_tests : remove_objects $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_BASIC_TESTS)      -lcunit -o $(BIN)/run_basic_tests_$(MAXK)$(READ_PAIR_SUFFIX) $(BASIC_TESTS_OBJ)

run_hash_table_tests : remove_objects $(HASH_TABLE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_HASH_TABLE_TESTS) -lcunit -o $(BIN)/run_hash_table_tests_$(MAXK) $(HASH_TABLE_TESTS_OBJ) 

run_graph_tests : remove_objects $(GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(LINKOPT) $(CFLAGS_CUNIT)  -o $(BIN)/run_graph_tests_$(MAXK) $(GRAPH_TESTS_OBJ) -lcunit 



tests: remove_objects run_basic_tests run_hash_table_tests run_graph_tests

# Cleaning rules
.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf obj

remove_objects:
	rm -rf obj

# Pattern rules

# cortex_con
obj/cortex_con/%.o : src/cortex_con/%.c include/cortex_con/%.h
	mkdir -p obj/cortex_con;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/cortex_con/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/cortex_con;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/cortex_con/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_value.h
	mkdir -p obj/cortex_con;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/cortex_con/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p obj/cortex_con;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $< -o $@

obj/cortex_con/%.o : src/cortex_con/%.c 
	mkdir -p obj/cortex_con;  $(CC) $(CFLAGS_CORTEX_CON) $(OPT) -c $? -o $@

# Unit tests
obj/test/%.o : src/test/cortex_con/%.c include/test/cortex_con/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_CORTEX_CON_TESTS) $(OPT) -c $< -o $@	

obj/test/run_dB_graph_tests.o : src/test/cortex_con/run_dB_graph_tests.c
	mkdir -p obj/test/; $(CC) $(CFLAGS_CORTEX_CON_TESTS) $(OPT) -c $< -o $@	

obj/test/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/basic/%.c include/test/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/hash_table/%.c include/test/hash_table/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

# Utils
obj/util/bubbleparse.o : src/util/bubbleparse.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_BUBBLEPARSE) $(OPT) -c $? -o $@

obj/util/graphout.o : src/util/graphout.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@

obj/util/filter_reads.o : src/util/filter_reads.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_FILTERREADS) $(OPT) -c $? -o $@

obj/util/graph_formats.o : src/util/graph_formats.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@

obj/util/graph_tools.o : src/util/graph_tools.c
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@

obj/util/node_queue.o : src/util/node_queue.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@

obj/util/kmerinfo.o : src/util/kmerinfo.c 
	mkdir -p obj/util;  $(CC) $(CFLAGS_GRAPHOUT) $(OPT) -c $? -o $@

obj/util/reformat_log.o : src/util/reformat_log.c
	mkdir -p obj/util;  $(CC) $(CFLAGS_LOGPARSE) $(OPT) -c $? -o $@

obj/util/cortex_log_parse.o : src/util/cortex_log_parse.c
	mkdir -p obj/util;  $(CC) $(CFLAGS_LOGPARSE) $(OPT) -c $? -o $@

obj/util/solid_converter.o : src/util/solid_converter.c
	mkdir -p obj/util;  $(CC) $(CFLAGS_LOGPARSE) $(OPT) -c $? -o $@
	
obj/util/subsampler.o : src/util/subsampler.c 
		mkdir -p obj/util;  $(CC) $(CFLAGS_BUBBLEPARSE) $(OPT) -c $? -o $@

# Stats
obj/stats/%.o : src/stats/%.c include/stats/%.h
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_STATS) $(OPT) -c $< -o $@

obj/stats/kmer_stats.o : src/stats/kmer_stats.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_STATS) $(OPT) -c $< -o $@

obj/stats/demultiplexer.o : src/stats/demultiplexer.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_DEMULTIPLEXER) $(OPT) -c $< -o $@

obj/stats/kmer_hash_build.o : src/stats/kmer_hash_build.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_HASH_BUILD) $(OPT) -c $< -o $@
	
obj/stats/kmer_filter.o : src/stats/kmer_filter.c
	mkdir -p obj/stats;  $(CC) $(CFLAGS_KMER_HASH_BUILD) $(OPT) -c $< -o $@


libctx.dylib:$(LIBRARY_OBJ)
	ld  -dylib -dynamic $(LIBRARY_OBJ) -o $(LIB)/libctx_$(BIN_SUFFIX).dylib -lc -lz ; 

libctx.so.1:$(LIBRARY_OBJ)
	$(CC) -shared -Wl,-soname,libctx.so -o $(LIB)/libctx_$(BIN_SUFFIX).so.1 $(LIBRARY_OBJ) -lc -lz

dylib:
	mkdir -p $(BIN);
	@$(MAKE) clean; \
	case `uname` in \
		Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libctx.so.1;; \
		Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libctx.dylib;; \
			*) echo 'Unknown OS';; \
	esac

