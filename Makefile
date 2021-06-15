UNAME_S := $(shell uname -s)

SRC = src/colord

ifeq ($(UNAME_S),Darwin)
	CC=/usr/local/bin/g++-10
	CFLAGS = -Wall -O3 -std=c++17 -static-libgcc -static-libstdc++ -pthread
	CLINK = -Wall -O3 -std=c++17 -static-libgcc -static-libstdc++ -lpthread

	CFLAGS_KMC = -Wall -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11
	CLINK_KMC = -lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11

else
	CC=g++
	CFLAGS = -Wall -O3 -std=c++17 -static -Wl,--whole-archive -lstdc++fs -lpthread -Wl,--no-whole-archive
	CLINK = -Wall -O3 -std=c++17 -static -Wl,--whole-archive -lstdc++fs -lpthread -Wl,--no-whole-archive	

	CFLAGS_KMC = -Wall -O3 -m64 -static-libgcc -static-libstdc++ -fopenmp -pthread -std=c++11
	CLINK_KMC = -lm -fopenmp -static-libgcc -static-libstdc++ -O3 -pthread -std=c++11
endif





BIN_DIR = bin

OBJS = \
$(SRC)/utils.o \
$(SRC)/timer.o \
$(SRC)/stats_collector.o \
$(SRC)/reads_sim_graph.o \
$(SRC)/quality_coder.o \
$(SRC)/quality_coder_impl.o \
$(SRC)/pooled_threads.o \
$(SRC)/main.o \
$(SRC)/in_reads.o \
$(SRC)/id_coder.o \
$(SRC)/filter_kmers.o \
$(SRC)/entr_header.o \
$(SRC)/encoder.o \
$(SRC)/dna_coder.o \
$(SRC)/decompression.o \
$(SRC)/info.o \
$(SRC)/count_kmers.o \
$(SRC)/compression.o \
$(SRC)/basic_coder.o \
$(SRC)/arg_parse.o \
$(SRC)/archive.o \
$(SRC)/reference_genome.o \
$(SRC)/libs/edlib/edlib.o \
$(SRC)/libs/kmc_api/kmc_file.o \
$(SRC)/libs/kmc_api/kmer_api.o \
$(SRC)/libs/kmc_api/mmer.o 

COBJS = \
$(SRC)/libs/md5/md5.o

ifeq ($(UNAME_S),Darwin)
    LIBS = \
	$(SRC)/../common/libs/zlib/libz.mac.a \
	$(SRC)/libs/mimalloc/libmimalloc.mac.a \
	$(SRC)/libs/count_kmers/libfiltering_kmc.mac.a \
	$(SRC)/libs/count_kmers/libbz2.1.0.5.dylib

	LIB_FILTERING_KMC = $(SRC)/libs/count_kmers/libfiltering_kmc.mac.a
else
	LIBS = \
	$(SRC)/../common/libs/zlib/libz.a \
	$(SRC)/libs/mimalloc/libmimalloc.a \
	$(SRC)/libs/count_kmers/libfiltering_kmc.a \
	$(SRC)/libs/count_kmers/libbz2.a

	LIB_FILTERING_KMC = $(SRC)/libs/count_kmers/libfiltering_kmc.a
endif


$(BIN_DIR)/colord: $(OBJS) $(COBJS) $(LIB_FILTERING_KMC)
	-mkdir -p $(BIN_DIR)
	$(CC) $(CLINK) -o $@ $^ $(LIBS)

$(OBJS): %.o: %.cpp
	$(CC) $(CFLAGS)  -I $(SRC)/../common/libs/zlib -I $(SRC)/libs/kmc_api -I $(SRC)/libs/edlib -I $(SRC)/libs/CLI11 -c $< -o $@

$(COBJS): %.o: %.c
	$(CC) $(CFLAGS)  -I $(SRC)/../common/libs/zlib -I $(SRC)/libs/kmc_api -I $(SRC)/libs/edlib -I $(SRC)/libs/CLI11 -c $< -o $@


KMC_MAIN_DIR = src/filtering-KMC

KMC_OBJS = \
$(KMC_MAIN_DIR)/kmer_counter.o \
$(KMC_MAIN_DIR)/mmer.o \
$(KMC_MAIN_DIR)/mem_disk_file.o \
$(KMC_MAIN_DIR)/rev_byte.o \
$(KMC_MAIN_DIR)/bkb_writer.o \
$(KMC_MAIN_DIR)/cpu_info.o \
$(KMC_MAIN_DIR)/bkb_reader.o \
$(KMC_MAIN_DIR)/fastq_reader.o \
$(KMC_MAIN_DIR)/timer.o \
$(KMC_MAIN_DIR)/develop.o \
$(KMC_MAIN_DIR)/kb_completer.o \
$(KMC_MAIN_DIR)/kb_storer.o \
$(KMC_MAIN_DIR)/kmer.o \
$(KMC_MAIN_DIR)/splitter.o \
$(KMC_MAIN_DIR)/kb_collector.o

ifeq ($(UNAME_S),Darwin)
	RADULS_OBJS = 
else
	RADULS_OBJS = \
	$(KMC_MAIN_DIR)/raduls_sse2.o \
	$(KMC_MAIN_DIR)/raduls_sse41.o \
	$(KMC_MAIN_DIR)/raduls_avx2.o \
	$(KMC_MAIN_DIR)/raduls_avx.o

endif

$(KMC_OBJS): %.o: %.cpp
	$(CC) $(CFLAGS_KMC) -I $(SRC)/../common/libs/zlib -c $< -o $@




$(KMC_MAIN_DIR)/raduls_sse2.o: $(KMC_MAIN_DIR)/raduls_sse2.cpp
	$(CC) $(CFLAGS_KMC) -msse2 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_sse41.o: $(KMC_MAIN_DIR)/raduls_sse41.cpp
	$(CC) $(CFLAGS_KMC) -msse4.1 -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx.o: $(KMC_MAIN_DIR)/raduls_avx.cpp
	$(CC) $(CFLAGS_KMC) -mavx -c $< -o $@
$(KMC_MAIN_DIR)/raduls_avx2.o: $(KMC_MAIN_DIR)/raduls_avx2.cpp
	$(CC) $(CFLAGS_KMC) -mavx2 -c $< -o $@
	


$(LIB_FILTERING_KMC): $(KMC_OBJS) $(RADULS_OBJS)
	-mkdir -p $(BIN_DIR)
	ar rcs $@ $^

clean:
	-rm -f $(SRC)/*.o
	-rm -f $(SRC)/libs/edlib/edlib.o
	-rm -f $(SRC)/libs/kmc_api/kmc_file.o
	-rm -f $(SRC)/libs/kmc_api/kmer_api.o
	-rm -f $(SRC)/libs/kmc_api/mmer.o
	-rm -f $(KMC_MAIN_DIR)/*.o
	-rm -f $(LIB_FILTERING_KMC)
	-rm -rf bin	


all:	$(BIN_DIR)/colord
