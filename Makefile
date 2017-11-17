all		: metacortex

ifndef CC
  CC = gcc
endif

ifdef MAC
    $(warning On MacOS, make sure you set CC to point to your GCC binary - LLVM won't compile MetaCortex)
    $(warning     e.g. export CC=/usr/local/Cellar/gcc/7.1.0/bin/gcc-7)
endif

BIN = bin
LIB = lib

ifndef MAXK
 MAXK = 31
endif

ifeq ($(shell expr $(MAXK) \< 32),1)
BITFIELDS = 1
BIN_SUFFIX=31
else ifeq ($(shell expr $(MAXK) \< 64),1)
BITFIELDS = 2
BIN_SUFFIX=63
else ifeq ($(shell expr $(MAXK) \< 96),1)
BITFIELDS = 3
BIN_SUFFIX=95
else ifeq ($(shell expr $(MAXK) \< 128),1)
BITFIELDS = 4
BIN_SUFFIX=127
else ifeq ($(shell expr $(MAXK) \< 161),1)
BITFIELDS = 5
BIN_SUFFIX=159
else ifeq ($(shell expr $(MAXK) \< 193),1)
BITFIELDS = 6
BIN_SUFFIX=191
else ifeq ($(shell expr $(MAXK) \< 224),1)
BITFIELDS = 7
BIN_SUFFIX=223
else ifeq ($(shell expr $(MAXK) \< 256),1)
BITFIELDS = 8
BIN_SUFFIX=255
else
BITFIELDS = 1
BIN_SUFFIX=31
endif

#BIN_SUFFIX = $(MAXK)

# Main program includes
IDIR_BASIC =include/basic
IDIR_UTIL =include/util
IDIR_HASH  =include/hash_table
IDIR_CORTEX = include/cortex
IDIR_ALIGNMENT = include/alignment
IDIR_STATS = include/stats

# Test code includes
IDIR_BASIC_TESTS =include/test/basic
IDIR_HASH_TABLE_TESTS =include/test/hash_table
IDIR_CORTEX_TESTS=include/test/cortex

#Default CUnit installation path in some Linux distributions.
IDIR_CUNIT=/usr/local/include/CUnit/

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
    MAC=1
endif

# Mac OS X specific. Assuming CUnit is installed with MacPorts
ifdef MAC
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
OPT		= $(ARCH) -Wall -O3 -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g
#-Wno-duplicate-decl-specifier

ifdef DEBUG
OPT	= $(ARCH) -Wall -O0 -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -g -pthread
endif

ifdef DEBUG_PRINT_LABELS
 OPT += -DDEBUG_PRINT_LABELS
endif

metacortex: OPT += -DNUMBER_OF_COLOURS=1
graphout:	OPT += -DNUMBER_OF_COLOURS=2
kmerinfo:	OPT += -DNUMBER_OF_COLOURS=2

# Include dirs
CFLAGS_METACORTEX = -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX) -I$(IDIR_UTIL)
CFLAGS_GRAPHOUT = -I$(IDIR_BASIC) -I$(IDIR_CORTEX) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_FILTERREADS = -I$(IDIR_BASIC) -I$(IDIR_CORTEX) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_KMERINFO = -I$(IDIR_BASIC) -I$(IDIR_CORTEX) -I$(IDIR_HASH)
CFLAGS_LOGPARSE = -I$(IDIR_BASIC) -I$(IDIR_CORTEX) -I$(IDIR_HASH) -I$(IDIR_UTIL)
CFLAGS_HASH_TABLE_TESTS = -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX) -I$(IDIR_HASH_TABLE_TESTS)
CFLAGS_METACORTEX_TESTS	= -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX) -I$(IDIR_CORTEX_TESTS)
CFLAGS_BASIC_TESTS = -I$(IDIR_CUNIT) $(CFLAGS_CUNIT) -I$(IDIR_BASIC) -I$(IDIR_BASIC_TESTS)

# Program objects
METACORTEX_OBJ = obj/cortex/file_format.o obj/cortex/flags.o obj/cortex/cleaning.o obj/cortex/path.o obj/cortex/perfect_path.o obj/cortex/branches.o obj/cortex/y_walk.o obj/cortex/cmd_line.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_value.o obj/cortex/hash_table.o obj/cortex/dB_graph.o obj/cortex/file_reader.o obj/cortex/metacortex.o obj/cortex/logger.o obj/cortex/metagraphs.o obj/cortex/coverage_walk.o obj/util/node_queue.o obj/cortex/graph_stats.o obj/cortex/bubble_find.o obj/cortex/report_output.o
KMERINFO_OBJ =  obj/cortex/flags.o obj/cortex/path.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_table.o obj/cortex/file_reader.o obj/util/kmerinfo.o obj/cortex/logger.o obj/cortex/hash_value.o
GRAPHOUT_OBJ = obj/cortex/flags.o obj/cortex/path.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_table.o obj/cortex/file_reader.o obj/cortex/dB_graph.o obj/util/graphout.o obj/cortex/perfect_path.o obj/cortex/logger.o obj/cortex/hash_value.o obj/util/graph_formats.o obj/util/node_queue.o obj/cortex/cleaning.o obj/cortex/coverage_walk.o obj/util/graph_tools.o obj/cortex/file_format.o
FILTERREADS_OBJ = obj/util/filter_reads.o obj/cortex/flags.o obj/cortex/path.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_table.o obj/cortex/file_reader.o obj/cortex/dB_graph.o obj/cortex/perfect_path.o obj/cortex/logger.o obj/cortex/hash_value.o obj/util/node_queue.o obj/cortex/cleaning.o obj/cortex/file_format.o
BASIC_TESTS_OBJ	= obj/test/binary_kmer.o obj/test/seq.o obj/test/test_binary_kmer.o obj/test/test_seq.o obj/test/run_basic_tests.o  obj/cortex/logger.o
HASH_TABLE_TESTS_OBJ = obj/cortex/flags.o obj/test/run_hash_table_tests.o obj/cortex/element.o obj/cortex/hash_value.o obj/cortex/hash_table.o obj/test/test_hash.o obj/cortex/binary_kmer.o  obj/cortex/seq.o obj/cortex/logger.o
GRAPH_TESTS_OBJ = obj/cortex/branches.o obj/cortex/file_format.o obj/test/test_dB_graph.o obj/cortex/logger.o  obj/cortex/cleaning.o  obj/cortex/perfect_path.o obj/cortex/path.o obj/cortex/flags.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_value.o obj/cortex/hash_table.o obj/cortex/dB_graph.o obj/cortex/file_reader.o obj/cortex/y_walk.o  obj/test/test_file_reader.o obj/test/test_graph_element.o obj/test/run_dB_graph_tests.o

#Library objects
LIBRARY_OBJ =  obj/cortex/file_format.o obj/cortex/analysis.o obj/cortex/flags.o obj/cortex/cleaning.o obj/cortex/path.o obj/cortex/perfect_path.o obj/cortex/branches.o obj/cortex/y_walk.o obj/cortex/cmd_line.o obj/cortex/binary_kmer.o obj/cortex/seq.o obj/cortex/element.o obj/cortex/hash_value.o obj/cortex/hash_table.o obj/cortex/dB_graph.o obj/cortex/file_reader.o obj/cortex/metacortex.o obj/cortex/logger.o obj/cortex/metagraphs.o obj/cortex/coverage_walk.o obj/util/node_queue.o obj/cortex/graph_stats.o obj/cortex/bubble_find.o obj/cortex/report_output.o

# Main rules
metacortex : remove_objects $(METACORTEX_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(OPT_COLS) -o $(BIN)/metacortex_k$(BIN_SUFFIX) $(METACORTEX_OBJ) -lm

kmerinfo: remove_objects $(KMERINFO_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -o $(BIN)/kmerinfo $(KMERINFO_OBJ) -lm

graphout:remove_objects $(GRAPHOUT_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -o $(BIN)/graphout_$(BIN_SUFFIX) $(GRAPHOUT_OBJ) -lm

filterreads:remove_objects $(FILTERREADS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -o $(BIN)/filterreads_$(BIN_SUFFIX) $(FILTERREADS_OBJ) -lm

run_basic_tests : remove_objects $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_BASIC_TESTS)      -lcunit -o $(BIN)/run_basic_tests_$(BIN_SUFFIX)$(READ_PAIR_SUFFIX) $(BASIC_TESTS_OBJ)

run_hash_table_tests : remove_objects $(HASH_TABLE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_HASH_TABLE_TESTS) -lcunit -o $(BIN)/run_hash_table_tests_$(BIN_SUFFIX) $(HASH_TABLE_TESTS_OBJ)

run_graph_tests : remove_objects $(GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(LINKOPT) $(CFLAGS_CUNIT)  -o $(BIN)/run_graph_tests_$(BIN_SUFFIX) $(GRAPH_TESTS_OBJ) -lcunit

tests: remove_objects run_basic_tests run_hash_table_tests run_graph_tests

# Cleaning rules
.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf obj

remove_objects:
	rm -rf obj

# Pattern rules

# Metacortex
obj/cortex/%.o : src/cortex/%.c include/cortex/%.h
	mkdir -p obj/cortex;  $(CC) $(CFLAGS_METACORTEX) $(OPT) -c $< -o $@

obj/cortex/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/cortex;  $(CC) $(CFLAGS_METACORTEX) $(OPT) -c $< -o $@

obj/cortex/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_value.h
	mkdir -p obj/cortex;  $(CC) $(CFLAGS_METACORTEX) $(OPT) -c $< -o $@

obj/cortex/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p obj/cortex;  $(CC) $(CFLAGS_METACORTEX) $(OPT) -c $< -o $@

obj/cortex/%.o : src/cortex/%.c
	mkdir -p obj/cortex;  $(CC) $(CFLAGS_METACORTEX) $(OPT) -c $? -o $@

# Unit tests
obj/test/%.o : src/test/cortex/%.c include/test/cortex/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_METACORTEX_TESTS) $(OPT) -c $< -o $@

obj/test/run_dB_graph_tests.o : src/test/cortex/run_dB_graph_tests.c
	mkdir -p obj/test/; $(CC) $(CFLAGS_METACORTEX_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/basic/%.c include/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/basic/%.c include/test/basic/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/test/%.o : src/test/hash_table/%.c include/test/hash_table/%.h
	mkdir -p obj/test/; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

# Utils
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
