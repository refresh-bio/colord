# CoLoRd - Compressing long reads

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/CoLoRd/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/CoLoRd/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/colord.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/colord)
[![C/C++ CI](https://github.com/refresh-bio/CoLoRd-dev/workflows/C/C++%20CI/badge.svg)](https://github.com/refresh-bio/CoLoRd-dev/actions)
[![License](https://anaconda.org/bioconda/colord/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

A versatile compressor of third generation sequencing reads.

## Quick start

```bash
git clone https://github.com/refresh-bio/CoLoRd
cd colord && make
cd bin

INPUT=./../test

# default compression presets (lossy quality, memory priority)
./colord compress-ont ${INPUT}/M.bovis.fastq ont.default 			# Oxford Nanopore
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

```

## Installation and configuration

CoLoRd comes with a set of [precompiled binaries](https://github.com/refresh-bio/CoLoRd/releases) for Windows, Linux, and OS X. They can be found under Releases tab. 
The software is also available on [Bioconda](https://anaconda.org/bioconda/colord):
```
conda install -c bioconda colord
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda).
CoLoRd can be also built from the sources distributed as:

* Visual Studio 2019 solution for Windows,
* MAKE project (G++ 8.4 required) for Linux and macOS.

To install G++ under under macOS, one can use *Homebrew* package manager:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install gcc
```
Before running CoLoRd on macOS, the current limit of file descriptors should be increased:
```
ulimit -n 2048
```

## Usage

### Compression

`colord <mode> [options] <input> <archive>`

Modes:
* `compress-ont` - compress Oxford Nanopore reads,
* `compress-pbhifi` - compress PacBio HiFi reads,
* `compress-pbraw` - compress PacBio CLR/subreads.

Positionals: 
* `input` - input FASTQ/FASTA path (gzipped or not),
* `output` - archive path. 

Options:
* `-h, --help` - print help
* `-k, --kmer-len` - *k*-mer length, (15-28, default: auto adjust)
* `-t, --threads` - number of threads (default: 12)
* `-p, --priority` - compression priority:  `memory`, `balanced`, `ratio` (default: `memory`)
* `-q, --qual` - quality compression mode: 
	* `org` - original,
	* `none` - discard (Q0 for all bases),
	* `avg` - average over entire file,
	* `2-fix`,`4-fix`,`5-fix` - 2/4/5 bins with fixed representatives,
	* `2-avg`,`4-avg`,`5-avg` - 2/4/5 bins with averages as representatives; default value depends on the mode (`4-avg` for `ont`, `5-avg` for `pbhifi`, `none` for `pbraw`),                           
* `-T, --qual-thresholds` - quality thresholds:
	* single value for `2-fix`/`2-avg` (default: 7),
	* three values for `4-fix`/`4-avg` (default: 7 14 26),
	* four values for `4-fix`/`4-avg` (default: 7 14 26 93),
	* not allowed for `avg`, `org` and `none` modes,
* `-D, --qual-values` - bin representatives for decompression,
   * single value for `none` mode (default: 0),
   * two values for `2-fix` mode (default: 1 13),
   * four values for `4-fix` mode (default: 3 10 18 35),
   * five values for `5-fix` mode (default: 3 10 18 35 93),
   * not allowed for `avg`, `org`, `2-avg`, `4-avg` and `5-avg` modes,
* `-G, --reference-genome` - optional reference genome path (multi-FASTA gzipped or not), it enables reference-based mode which provides better compression ratios,
* `-s, --store-reference` - stores the reference genome in the archive, use only with `-G` flag,
* `-v, --verbose` - verbose mode.

Advanced options (default values may depend on the mode - please run `colord --help <mode>` to get the details):                             
* `-a, --anchor-len` - anchor len (default: auto adjust),
* `-L, --Lowest-count` - minimal *k*-mer count,
* `-H, --Highest-count` - maximal *k*-mer count,
* `-f, --filter-modulo` - k-mers for which *hash(k-mer) mod f != 0* will be filtered out before graph building,
* `-c, --max-candidates` - maximal number of reference reads considered as reference,
* `-e, --edit-script-mult` - multipier for predicted cost of storing read part as edit script,
* `-r, --max-recurence-level` - maximal level of recurence when considering alternative reference reads,
* `--min-to-alt` - minimum length of encoding part to consider using alternative read,
* `--min-mmer-frac` - if *A* is set of m-mers in encode read R then read is refused from encoding if *|A| < min-mmer-frac * len(R)*,
* `--min-mmer-force-enc` - if *A* is set of m-mers in encode read R then read is accepted to encoding always if *|A| > min-mmer-force-enc * len(R)*,
* `--max-matches-mult` - if the number of matches between encode read *R* and reference read is *r*, then read is refused from encoding if *r > max-matches-mult * len(R)*,
* `--fill-factor-filtered-kmers` - fill factor of filtered *k*-mers hash table,
* `--fill-factor-kmers-to-reads` - fill factor of *k*-mers to reads hash table,
* `--min-anchors` - if number of anchors common to encode read and reference candidate is lower than minAnchors candidate is refused,
* `-i, --identifier` header compression mode - `main`/`none`/`org` (default: `org`),                        
* `-R, --Ref-reads-mode` - reference reads mode: `all`/`sparse` (default: `sparse`),                             
* `-g, --sparse-range` - sparse mode range. The propability of reference read acceptance is *1 / pow(id/range_reads, exponent)*, where range_reads is determined based on the number of symbols, which in turn is determined by the number of trusted unique *k*-mers (estimated genome length) multiplied by the value of this parameter,
* `-x, --sparse-exponent` - sparse mode exponent.

### Decompression

`colord decompress [options] <archive> <output>`

Positionals:
* `input` - archive path,
* `output` - output file path.

Options:
* `-h, --help` - print help,
* `-G, --reference-genome` - optional reference genome path (multi-FASTA gzipped or not), required for reference-based archives with no reference genome embedded (`-G` compression without `-s` switch),
* `-v, --verbose` - verbose mode.


### Archive information

`colord info <archive>`


## Citing
[paper](link)




