# !/bin/bash

cd bin
INPUT=./../test

# default compression presets (lossy quality, memory priority)
./colord compress-ont ${INPUT}/M.bovis.fastq ont.default 		# Oxford Nanopore
./colord compress-pbhifi ${INPUT}/D.melanogaster.fastq hifi.default	# PacBio HiFi 
./colord compress-pbraw ${INPUT}/A.thaliana.fastq clr.default 		# PacBio CLR/subreads

# print ONT archive information and decompress
./colord info ont.default
./colord decompress ont.default ont.fastq

# compress HiFi reads preserving original quality levels
./colord compress-pbhifi -q org ${INPUT}/D.melanogaster.fastq hifi.lossless

# compress CLR reads with ratio priority using 48 threads
./colord compress-pbraw -p ratio -t 48 ${INPUT}/A.thaliana.fastq clr.ratio

# compress ONT reads w.r.t. reference genome (embed the reference in the archive)
./colord compress-ont -G ${INPUT}/M.bovis-reference.fna -s ${INPUT}/M.bovis.fastq ont.refbased

# decompress the reference-based archive
./colord decompress ont.refbased ont.refbased.fastq
