# CoLoRd - Compressing long reads

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/colord/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/colord/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/colord.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/colord)
[![GitHub Actions CI](../../actions/workflows/main.yml/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A versatile compressor of third generation sequencing reads.

## Quick start

```bash
git clone https://github.com/refresh-bio/colord
cd colord && make
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

```

## Installation and configuration

CoLoRd comes with a set of [precompiled binaries](https://github.com/refresh-bio/colord/releases) for Windows, Linux, and OS X. They can be found under Releases tab. 
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
brew install gcc@10
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

#### Hints
While the number of CoLoRd parameters is large, in most cases the default values will work just fine.
In terms of compression, there is always a trade off between compression ratio and resource requirements (mainly memory and compute time).
If the default behavior of CoLoRd is insufficient, the first attempt should be the change of compression priority mode (```-p``` parameter). 
The compression priority modes aggregate multiple other parameters influencing compression ratio.
There are the following priority modes (ordered increasingly w.r.t. the compression efficiency and resource requirements): 

 * ```memory``` 
 * ```balanced``` 
 * ```ratio``` 

The ```memory``` priority mode is the default.

Quality scores have a high impact on the compression. They are hard to compress due to their nature and, at the same time (as presented in the paper) their resolution can be safely reduced without affecting downstream analyses. For this reason, in each  priority mode, the quality scores are compressed lossy. If it is required to keep the original quality scores, one should use ```-q org```. Note, that there exist several other quality compression modes (see the paper).

Here are compression results for a large set of human reads [NA12878](http://s3.amazonaws.com/nanopore-human-wgs/rel6/rel_6.fastq.gz) with a total size of 268,305,314,354 bytes.

|                                            | Lossy           | Lossless        |
| ------------------------------------------ | --------------- | --------------- |
| Compressed in ```memory``` mode size [B]   | 42,120,596,486  | 105,807,350,384 |
| Compressed in ```balanced``` mode size [B] | 39,833,878,505  | 103,367,993,362 |
| Compressed in ```ratio``` mode size [B]    | 38,832,714,102  | 101,305,368,675 |
| Time in ```memory``` mode [h:mm:ss]        | 1:12:42         | 1:26:02         |
| Time in ```balanced``` mode [h:mm:ss]      | 1:33:18         | 2:11:21         |
| Time in ```ratio``` mode [h:mm:ss]         | 3:18:46         | 4:57:09         |
| Memory in ```memory``` mode [KB]           | 13,715,168      | 14,341,128      |
| Memory in ```balanced``` mode [KB]         | 26,728,108      | 27,293,824      |
| Memory in ```ratio``` mode [KB]            | 97,922,208      | 99,133,548      |


If one wants to check how much CoLoRd can squeeze the input data regardless of the resource requirements, the ```ratio``` mode should be used.
If more control over execution is in demand, the remaining parameters may be configured. 
The simplest way to settle the direction without the need to understand the meaning of parameters is to display the defaults for a given compression priority mode with ```--help``` switch.
For example, let's say you want to find out if you should increase or decrease the ```-f``` parameter to improve the compression ratio while compressing ONT data.
You may run CoLoRd twice with the following parameters:
```
./colord compress-ont --help -p balanced
./colord compress-ont --help -p ratio
``` 
You will notice the default for ```-f``` is higher for ```balanced``` mode, which means lowering it will increase the compression ratio. The same approach may be applied for other parameters (```-L```, ```-H```, ```-c```, ```-r```, ```--min-to-alt```, etc.).

In the ```ratio``` priority mode all the input reads may serve as a reference to encode other reads. This will increase RAM usage, especially for large datasets. In the remaining modes, only part of the reads may serve as a reference. If needed ```-g``` and ```-x``` may be used.

The values for ```-k``` and ```-a``` parameters are auto-adjusted based on the size of the data to be compressed. The general rule is, the larger the input size is, the values of these parameters should be higher.



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

## API

CoLoRd comes with a C++ API allowing straightforward access to the existing archive. Below one can find an example of using API in the code.

```c++
#include "colord_api.h"
#include <iostream>

int main(int argc, char** argv) {
	try {
		colord::DecompressionStream stream("archive.colord");	// load a CoLoRd archive
		auto info = stream.GetInfo();				// get and print archive information
		std::cerr << "Archive info:\n\n";			//
		info.ToOstream(std::cerr);				//	

     		// iterate over records in the archive
		while (auto x = stream.NextRecord()) {
			if (info.isFastq) {
				std::cout << "@" << x.ReadHeader() << "\n";
				std::cout << x.Read() << "\n";
				std::cout << "+" << x.QualHeader() << "\n";
				std::cout << x.Qual() << "\n";
			} else {
				std::cout << ">" << x.ReadHeader() << "\n";
				std::cout << x.Read() << "\n";
			}
		}
	}
	catch (const std::exception& ex) {
		std::cerr << "Error: " << ex.what() << "\n";
		return -1;
	}	
	return 0;
}
```

### Compiling own code utilizing colord API
To use an API one needs to include ```colord_api.h``` header file and link against ```libcolord_api.a```. ```libcolord_api.a``` uses ```std::thread```s and zlib, so ```-lpthreads``` and ```-lz``` flags are needed for linking. For example, to compile and link the code above one could use the following command:
```
g++ -O3 $SRC_FILE -I$INCLUDE_DIR $LIB_DIR/libcolord_api.a -lz -lpthread -o example -no-pie
```
where
 * ```SRC_FILE``` is a path to a source code
 * ```INCLUDE_DIR``` is a path of the directory where ```colord_api.h``` file is (when one compiles ```colord``` from sources there is ```include``` directory created at the same location where ```Makefile``` is)
 * ```LIB_DIR``` is a path of the directory where ```libcolord_api.a``` file is (when one compiles ```colord``` from sources there is ```bin``` directory created at the same location where ```Makefile``` is, it contains (among others) ```libcolord_api.a```)


## Citing
[Marek Kokot, Adam Gudyś, Heng Li, Sebastian Deorowicz (2021) CoLoRd: Compressing long reads. *bioRxiv* 2021.07.17.452767; doi: https://doi.org/10.1101/2021.07.17.452767](https://doi.org/10.1101/2021.07.17.452767)




