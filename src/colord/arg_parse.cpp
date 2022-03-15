/*******************************************************************************
 
 CoLoRd 
 Copyright (C) 2021, M. Kokot, S. Deorowicz, and A. Gudys
 https://github.com/refresh-bio/CoLoRd

 This program is free software: you can redistribute it and/or modify it under 
 the terms of the GNU General Public License as published by the Free Software 
 Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY 
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this 
 program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#include "defs.h"
#include "arg_parse.h"
#include "utils.h"
//for some reason CLI11 causes link error (linux) or runtimer error (mac os) when filesystem is used
#ifdef __clang__
#define CLI11_HAS_FILESYSTEM 0
#endif
#include "CLI11.hpp"
#include "compression.h"
#include "decompression.h"
#include "info.h"
#include <string>

struct CDefaultQualBinThr
{
    static const inline std::vector<uint32_t> Fwd = { 7 };
    static const inline std::vector<uint32_t> Rev = { 1, 13 };
};

struct CDefaultQualQuadThr
{
    static const inline std::vector<uint32_t> Fwd = { 7, 14, 26 };
    static const inline std::vector<uint32_t> Rev = { 3, 10, 18, 35 };
};

struct CDefaultQualQuinaryThr
{
    static const inline std::vector<uint32_t> Fwd = { 7, 14, 26, 93 };
    static const inline std::vector<uint32_t> Rev = { 3, 10, 18, 35, 93 };
};

struct CDefaultQualNone
{
    static const inline std::vector<uint32_t> Fwd = { };
    static const inline std::vector<uint32_t> Rev = { 0 };
};

struct CDefaultQualAvg
{
    static const inline std::vector<uint32_t> Fwd = { };
    static const inline std::vector<uint32_t> Rev = { };
};

struct CDefaultQualOrg
{
    static const inline std::vector<uint32_t> Fwd = { };
    static const inline std::vector<uint32_t> Rev = { };
};

struct CDefaultQualBinAvg
{
    static const inline std::vector<uint32_t> Fwd = { 7 };
    static const inline std::vector<uint32_t> Rev = { };
};

struct CDefaultQualQuadAvg
{
    static const inline std::vector<uint32_t> Fwd = { 7, 14, 26 };
    static const inline std::vector<uint32_t> Rev = { };
};

struct CDefaultQualQuinaryAvg
{
    static const inline std::vector<uint32_t> Fwd = { 7, 14, 26, 93 };
    static const inline std::vector<uint32_t> Rev = { };
};

/********************************************************/
/*                    ONT defaults                      */
/********************************************************/
CCompressorParams compr_ONT_ratio_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::ONT;
    params.priority = CompressionPriority::Ratio;
    
    params.compressionLevel = 3;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 2;
    params.maxKmerCount = 120;
    params.filterHashModulo = 8;
    params.maxCandidates = 10;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 6;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    //    params.qualityComprMode = QualityComprMode::Original;
    params.qualityComprMode = QualityComprMode::QuadAverage;

    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::All;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}

CCompressorParams compr_ONT_balanced_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::ONT;
    params.priority = CompressionPriority::Balanced;

    params.compressionLevel = 2;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 3;
    params.maxKmerCount = 100;
    params.filterHashModulo = 9;
    params.maxCandidates = 8;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 5;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::QuadAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
    params.sparseMode_range_symbols = 2;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}


CCompressorParams compr_ONT_memory_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::ONT;
    params.priority = CompressionPriority::Memory;

    params.compressionLevel = 1;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 4;
    params.maxKmerCount = 80;
    params.filterHashModulo = 12;
    params.maxCandidates = 5;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 3;
    params.minPartLenToConsiderAltRead = 64;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);

    params.qualityComprMode = QualityComprMode::QuadAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}
/*------------------------------------------------------*/

/********************************************************/
/*                   PBRaw defaults                     */
/********************************************************/
CCompressorParams compr_PBRaw_ratio_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBRaw;
    params.priority = CompressionPriority::Ratio;

    params.compressionLevel = 3;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 2;
    params.maxKmerCount = 120;
    params.filterHashModulo = 8;
    params.maxCandidates = 10;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 6;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::None;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::All;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}

CCompressorParams compr_PBRaw_balanced_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBRaw;
    params.priority = CompressionPriority::Balanced;

    params.compressionLevel = 2;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 3;
    params.maxKmerCount = 100;
    params.filterHashModulo = 9;
    params.maxCandidates = 8;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 5;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::None;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
    params.sparseMode_range_symbols = 2;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}

CCompressorParams compr_PBRaw_memory_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBRaw;
    params.priority = CompressionPriority::Memory;

    params.compressionLevel = 1;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 4;
    params.maxKmerCount = 80;
    params.filterHashModulo = 12;
    params.maxCandidates = 5;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 3;
    params.minPartLenToConsiderAltRead = 64;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::None;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}
/*------------------------------------------------------*/

/********************************************************/
/*                  PBHiFi defaults                     */
/********************************************************/
CCompressorParams compr_PBHiFi_ratio_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBHiFi;
    params.priority = CompressionPriority::Ratio;

    params.compressionLevel = 3;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 2;
    params.maxKmerCount = 150;
    params.filterHashModulo = 20;
    params.maxCandidates = 12;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 6;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::QuinaryAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::All;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}

CCompressorParams compr_PBHiFi_balanced_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBHiFi;
    params.priority = CompressionPriority::Balanced;

    params.compressionLevel = 2;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 3;
    params.maxKmerCount = 120;
    params.filterHashModulo = 30;
    params.maxCandidates = 10;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 5;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::QuinaryAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
//    params.sparseMode_range_symbols = 2;
    params.sparseMode_range_symbols = 6;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
}

CCompressorParams compr_PBHiFi_memory_set_defaults()
{
    CCompressorParams params;
    params.dataSource = DataSource::PBHiFi;
    params.priority = CompressionPriority::Memory;

    params.compressionLevel = 2;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 3;
    params.maxKmerCount = 100;
    params.filterHashModulo = 40;
    params.maxCandidates = 8;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 5;
    params.minPartLenToConsiderAltRead = 48;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::QuinaryAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
//    params.sparseMode_range_symbols = 2;
    params.sparseMode_range_symbols = 3;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;
    
/*    CCompressorParams params;
    params.dataSource = DataSource::PBHiFi;
    params.priority = CompressionPriority::Memory;

    params.compressionLevel = 1;

    params.kmerLen = 0;
    params.anchorLen = 0;
    params.minKmerCount = 4;
    params.maxKmerCount = 80;
    params.filterHashModulo = 64;
    params.maxCandidates = 5;
    params.editScriptCostMultiplier = 1.0;
    params.maxRecurence = 3;
    params.minPartLenToConsiderAltRead = 64;
    params.minFractionOfMmersInEncode = 0.5;
    params.minFractionOfMmersInEncodeToAlwaysEncode = 0.9;
    params.maxMatchesMultiplier = 10;
    params.nThreads = std::max(std::thread::hardware_concurrency(), 1u);
    params.qualityComprMode = QualityComprMode::QuinaryAverage;
    params.headerComprMode = HeaderComprMode::Original;
    params.referenceReadsMode = ReferenceReadsMode::Sparse;
    params.sparseMode_range_symbols = 1;
    params.sparseMode_exponent = 1.0;
    params.minAnchors = 1;

    return params;*/
}
/*------------------------------------------------------*/



template<typename T>
void adjust_quality_mode_and_thresholds(CCompressorParams& params, const T& defaults, const std::string& mode_name)
{
    //set defaults if not set by params
    if (!params.qualityFwdThresholds.size())
        params.qualityFwdThresholds = defaults.Fwd;
    if (!params.qualityRevThresholds.size())
        params.qualityRevThresholds = defaults.Rev;

    if (params.qualityFwdThresholds.size() > defaults.Fwd.size())
    {
        std::cerr << "Warning: for '" + mode_name + "' quality compression mode expected number of quality thresholds is "
            + std::to_string(defaults.Fwd.size()) + ", but " + std::to_string(params.qualityFwdThresholds.size()) + " given. "
            "Last " + std::to_string(params.qualityFwdThresholds.size() - defaults.Fwd.size()) + " values will be ignored.\n";
        params.qualityFwdThresholds.resize(defaults.Fwd.size());
    }

    if (params.qualityRevThresholds.size() > defaults.Rev.size())
    {
        std::cerr << "Warning: for '" + mode_name + "' quality compression mode expected number of quality values is "
            + std::to_string(defaults.Rev.size()) + ", but " + std::to_string(params.qualityRevThresholds.size()) + " given. "
            "Last " + std::to_string(params.qualityRevThresholds.size() - defaults.Rev.size()) + " values will be ignored.\n";
        params.qualityRevThresholds.resize(defaults.Rev.size());
    }

    if (params.qualityFwdThresholds.size() < defaults.Fwd.size())
    {
        std::cerr << "Error:  for '" + mode_name + "' quality compression mode expected number of quality thresholds is "
            + std::to_string(defaults.Fwd.size()) + ", but " + std::to_string(params.qualityFwdThresholds.size()) + " given.\n";
        exit(1);
    }

    if (params.qualityRevThresholds.size() < defaults.Rev.size())
    {
        std::cerr << "Error: for '" + mode_name + "' quality compression mode expected number of quality values is "
            + std::to_string(defaults.Rev.size()) + ", but " + std::to_string(params.qualityRevThresholds.size()) + " given.\n";
        exit(1);
    }
}

void addPriorityParam(CLI::App& app, std::string& str)
{
    std::set<std::string> q_p{ "ratio", "balanced", "memory" };
    str = "memory"; //memory is default
    app.add_set("-p,--priority", str, q_p, "compression quality");
}

void set_compression_params(CLI::App* app, CInfo& info, const std::string& commandName, const std::string& commandDesc, CCompressorParams& comParams, bool hidden, std::vector<CLI::Option*>& toHideIfNoHelp)
{
    CLI::App* compParser = app->add_subcommand(commandName, commandDesc)->ignore_case();
    if (hidden)
        compParser->group("");

    //positionals    
    compParser->add_option("input", comParams.inputFilePath, "input FASTQ/FASTA path (gzipped or not)")->required(true)->check(CLI::ExistingFile);
    compParser->add_option("output", comParams.outputFilePath, "archive path")->required(true);

    // options
    toHideIfNoHelp.push_back(compParser->add_option("-k,--kmer-len", comParams.kmerLen, "k-mer length (default: auto adjust)")->check(CLI::Range(15, 28))); //TODO: maybe max should be lower (27, 28?)
    toHideIfNoHelp.push_back(compParser->add_option("-t,--threads", comParams.nThreads, "number of threads", true));
  
    addPriorityParam(*compParser, comParams.internal.priority);

    comParams.internal.qual_mode = qualityComprModeToString(comParams.qualityComprMode);

    std::set<std::string> q_p{ "org", "none", "2-fix", "4-fix", "5-fix", "avg", "2-avg", "4-avg", "5-avg" };
    toHideIfNoHelp.push_back(compParser->add_set("-q,--qual", comParams.internal.qual_mode, q_p, "quality compression mode \n"
        " * org - original, \n"
        " * none - discard (Q1 for all bases), \n"
        " * avg - average over entire file, \n"
        " * 2-fix, 4-fix, 5-fix - 2/4/5 bins with fixed representatives, \n"
        " * 2-avg, 4-avg, 5-avg - 2/4/5 bins with averages as representatives.",
        true));

    toHideIfNoHelp.push_back(compParser->add_option("-T,--qual-thresholds", comParams.qualityFwdThresholds,
        "quality thresholds, \n"
        " * single value for '2-fix' mode (default is" + vec_to_string(CDefaultQualBinThr::Fwd) + "), \n"
        " * single value for '2-avg' mode (default is" + vec_to_string(CDefaultQualBinAvg::Fwd) + "), \n"
        " * three values for '4-fix' mode (default is" + vec_to_string(CDefaultQualQuadThr::Fwd) + "), \n"
        " * three values for '4-avg' mode (default is" + vec_to_string(CDefaultQualQuadAvg::Fwd) + "), \n"
        " * four values for '5-fix' mode (default is" + vec_to_string(CDefaultQualQuinaryThr::Fwd) + "), \n"
        " * four values for '5-avg' mode (default is" + vec_to_string(CDefaultQualQuinaryAvg::Fwd) + "), \n"
        " * not allowed for 'avg', 'org' and 'none' modes", false));

    toHideIfNoHelp.push_back(compParser->add_option("-D,--qual-values", comParams.qualityRevThresholds,
        "\nquality values for decompression,\n"
        " * single value for 'none' mode (default is" + vec_to_string(CDefaultQualNone::Rev) + "), \n"
        " * two values for '2-fix' mode (default is" + vec_to_string(CDefaultQualBinThr::Rev) + "), \n"
        " * four values for '4-fix' mode (default is" + vec_to_string(CDefaultQualQuadThr::Rev) + "), \n"
        " * five values for '5-fix' mode (default is" + vec_to_string(CDefaultQualQuinaryThr::Rev) + "), \n"
        " * not allowed for 'avg', 'org', '2-avg', '4-avg' and '5-avg' modes", false));

    toHideIfNoHelp.push_back(compParser->add_option("-G,--reference-genome", comParams.refGenomePath, 
        "optional reference genome path (multi-FASTA gzipped or not), enables reference-based mode which provides better compression ratios")->check(CLI::ExistingFile));
    
    toHideIfNoHelp.push_back(compParser->add_flag("-s,--store-reference", comParams.storeRefGenome, 
        "stores the reference genome in the archive, use only with `-G` flag"));

    toHideIfNoHelp.push_back(compParser->add_flag("-v,--verbose", comParams.verbose, "verbose"));


    // pro options
    toHideIfNoHelp.push_back(compParser->add_option("-a,--anchor-len", comParams.anchorLen, "anchor len (default: auto adjust)"));
    toHideIfNoHelp.push_back(compParser->add_option("-L,--Lowest-count", comParams.minKmerCount, "minimal k-mer count", true));
    toHideIfNoHelp.push_back(compParser->add_option("-H,--Highest-count", comParams.maxKmerCount, "maximal k-mer count", true));
    toHideIfNoHelp.push_back(compParser->add_option("-f,--filter-modulo", comParams.filterHashModulo, "k-mers for which hash(k-mer) mod f != 0 will be filtered out before graph building", true));
    toHideIfNoHelp.push_back(compParser->add_option("-c,--max-candidates", comParams.maxCandidates, "maximal number of reference reads considered as reference", true)->check(CLI::PositiveNumber));
    toHideIfNoHelp.push_back(compParser->add_option("-e,--edit-script-mult", comParams.editScriptCostMultiplier, "multipier for predicted cost of storing read part as edit script", true)->check(CLI::PositiveNumber));
    toHideIfNoHelp.push_back(compParser->add_option("-r,--max-recurence-level", comParams.maxRecurence, "maximal level of recurence when considering alternative reference reads", true)->check(CLI::NonNegativeNumber));
    toHideIfNoHelp.push_back(compParser->add_option("--min-to-alt", comParams.minPartLenToConsiderAltRead, "minimum length of encoding part to consider using alternative read", true)->check(CLI::PositiveNumber));
    toHideIfNoHelp.push_back(compParser->add_option("--min-mmer-frac", comParams.minFractionOfMmersInEncode, "if A is set of m-mers in encode read R then read is refused from encoding if |A| < min-mmer-frac * len(R)", true)->check(CLI::NonNegativeNumber));
    toHideIfNoHelp.push_back(compParser->add_option("--min-mmer-force-enc", comParams.minFractionOfMmersInEncodeToAlwaysEncode, "if A is set of m-mers in encode read R then read is accepted to encoding always if |A| > min-mmer-force-enc * len(R)", true)->check(CLI::NonNegativeNumber));
    toHideIfNoHelp.push_back(compParser->add_option("--max-matches-mult", comParams.maxMatchesMultiplier, "if the number of matches between encode read R and reference read is r, then read is refused from encoding if r > max-matches-mult * len(R)", true)->check(CLI::NonNegativeNumber));
    toHideIfNoHelp.push_back(compParser->add_option("--fill-factor-filtered-kmers", comParams.fillFactorFilteredKmers, "fill factor of filtered k-mers hash table", true)->check(CLI::Range(0.1, 0.99)));
    toHideIfNoHelp.push_back(compParser->add_option("--fill-factor-kmers-to-reads", comParams.fillFactorKmersToReads, "fill factor of k-mers to reads hash table", true)->check(CLI::Range(0.1, 0.99)));
  
   
    toHideIfNoHelp.push_back(compParser->add_option("--min-anchors", comParams.minAnchors, "if number of anchors common to encode read and reference candidate is lower than minAnchors candidate is refused", true));

    comParams.internal.header_mode = headerComprModeToString(comParams.headerComprMode);
    
    std::set<std::string> h_p{ "org", "main", "none" };
    toHideIfNoHelp.push_back(compParser->add_set("-i,--identifier", comParams.internal.header_mode, h_p, "header compression mode", true));

    comParams.internal.reference_mode = referenceReadsModeToString(comParams.referenceReadsMode);
    
    std::set<std::string> ref_p{ "all", "sparse" };
    toHideIfNoHelp.push_back(compParser->add_set("-R,--Ref-reads-mode", comParams.internal.reference_mode, ref_p, "reference reads mode", true));

    toHideIfNoHelp.push_back(compParser->add_option("-g,--sparse-range", comParams.sparseMode_range_symbols, "sparse mode range. The propability of reference read acceptance is 1/pow(id/range_reads, exponent), where range_reads is determined based on the number of symbols, which in turn is determined by the number of trusted unique k-mers (estimated genome length) multiplied by the value of this parameter", true));
    toHideIfNoHelp.push_back(compParser->add_option("-x,--sparse-exponent", comParams.sparseMode_exponent, "sparse mode exponent", true));

    compParser->callback([&]() {

        comParams.priority = compressionPriorityFromString(comParams.internal.priority);

        comParams.headerComprMode = headerComprModeFromString(comParams.internal.header_mode);
        
        if (comParams.internal.qual_mode == "none")
        {
            comParams.qualityComprMode = QualityComprMode::None;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualNone{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "org")
        {
            comParams.qualityComprMode = QualityComprMode::Original;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualOrg{}, comParams.internal.qual_mode);
        }        
        else if (comParams.internal.qual_mode == "2-fix")
        {
            comParams.qualityComprMode = QualityComprMode::BinaryThreshold;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualBinThr{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "4-fix")
        {
            comParams.qualityComprMode = QualityComprMode::QuadThreshold;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualQuadThr{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "5-fix")
        {
            comParams.qualityComprMode = QualityComprMode::QuinaryThreshold;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualQuinaryThr{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "avg")
        {
            comParams.qualityComprMode = QualityComprMode::Average;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualAvg{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "2-avg")
        {
            comParams.qualityComprMode = QualityComprMode::BinaryAverage;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualBinAvg{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "4-avg")
        {
            comParams.qualityComprMode = QualityComprMode::QuadAverage;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualQuadAvg{}, comParams.internal.qual_mode);
        }
        else if (comParams.internal.qual_mode == "5-avg")
        {
            comParams.qualityComprMode = QualityComprMode::QuinaryAverage;
            adjust_quality_mode_and_thresholds(comParams, CDefaultQualQuinaryAvg{}, comParams.internal.qual_mode);
        }
        else
        {
            std::cerr << "Internal error\n";
            exit(1);
        }

        comParams.referenceReadsMode = referenceReadsModeFromString(comParams.internal.reference_mode);
        
        if (comParams.kmerLen && !comParams.anchorLen)
        {
            std::cerr << "Error: if -k,--kmer-len is set -a,--anchor-len also must be set\n";
            exit(1);
        }
        if (!comParams.kmerLen && comParams.anchorLen)
        {
            std::cerr << "Error: if -a,--anchor-len is set -k,--kmer-len also must be set\n";
            exit(1);
        }

        if (comParams.kmerLen && comParams.anchorLen)
        {
            if (comParams.anchorLen > comParams.kmerLen)
            {
                std::cerr << "Error: -a,--anchor-len must be less than or equal to -k,--kmer-len\n";
                exit(1);
            }
        }

        if (comParams.storeRefGenome && comParams.refGenomePath == "")
        {
            std::cerr << "Warning: -s has not effect if reference genome (-G) is not specified\n";
            comParams.storeRefGenome = false;
        }

        runCompression(comParams, info);
    });
}


std::string dataSourceToString(DataSource dataSource)
{
    switch (dataSource)
    {
    case DataSource::ONT:
        return "ONT";
    case DataSource::PBRaw:
        return "PBRaw";
    case DataSource::PBHiFi:
        return "PBHiFi";
    default:
        std::cerr << "Internal Error\n";
        exit(1);
    }
}

std::string qualityComprModeToString(QualityComprMode qualityComprMode)
{
    switch (qualityComprMode)
    {
    case QualityComprMode::None:
        return "none";
    case QualityComprMode::Original:
        return "org";
    case QualityComprMode::BinaryThreshold:
        return "2-fix";        
    case QualityComprMode::QuadThreshold:
        return "4-fix";        
    case QualityComprMode::QuinaryThreshold:
        return "5-fix";        
    case QualityComprMode::Average:
        return "avg";        
    case QualityComprMode::BinaryAverage:
        return "2-avg";        
    case QualityComprMode::QuadAverage:
        return "4-avg";        
    case QualityComprMode::QuinaryAverage:
        return "5-avg";        
    default:
        std::cerr << "Internal error\n";
        exit(1);
        break;
    }
}

HeaderComprMode headerComprModeFromString(const std::string& str)
{
    if (str == "org")
    {
        return HeaderComprMode::Original;
    }
    else if (str == "main")
    {
        return HeaderComprMode::Main;
    }
    else if (str == "none")
    {
        return HeaderComprMode::None;
    }
    else
    {
        std::cerr << "Internal error\n";
        exit(1);
    }    
}
std::string headerComprModeToString(HeaderComprMode headerComprMode)
{
    switch (headerComprMode)
    {
    case HeaderComprMode::Original:
        return "org";
        break;
    case HeaderComprMode::Main:
        return "main";
        break;
    case HeaderComprMode::None:
        return "none";
        break;
    default:
        std::cerr << "Internal error\n";
        exit(1);        
    }
}

CompressionPriority compressionPriorityFromString(const std::string& str)
{
    if (str == "ratio")
        return CompressionPriority::Ratio;
    else if (str == "balanced")
        return CompressionPriority::Balanced;
    else if (str == "memory")
        return CompressionPriority::Memory;
    else
    {
        std::cerr << "Internal error\n";
        exit(1);
    }
}
std::string compressionPriorityToString(CompressionPriority compressionPriority)
{
    switch (compressionPriority)
    {
    case CompressionPriority::Ratio:
        return "ratio";
    case CompressionPriority::Balanced:
        return "balanced";
    case CompressionPriority::Memory:
        return "memory";
    default:
        std::cerr << "Internal Error\n";
        exit(1);
    }
}

ReferenceReadsMode referenceReadsModeFromString(const std::string& str)
{
    if (str == "all")
        return ReferenceReadsMode::All;
    else if (str == "sparse")
        return ReferenceReadsMode::Sparse;
    else
    {
        std::cerr << "Internal error\n";
        exit(1);
    }
}
std::string referenceReadsModeToString(ReferenceReadsMode referenceReadsMode)
{
    switch (referenceReadsMode)
    {
    case ReferenceReadsMode::All:
        return "all";
        break;
    case ReferenceReadsMode::Sparse:
        return "sparse";
        break;
    default:
        std::cerr << "Internal error\n";
        exit(1);
        break;
    }
}

std::string vec_to_string(const std::vector<uint32_t>& vec)
{
    std::ostringstream oss;
    for (auto v : vec)
        oss << " " << v;
    return oss.str();
}


// fake parsing, only for checking if priority was set
CompressionPriority checkPriority(int argc, char** argv)
{
    CLI::App tmp{};
    tmp.set_help_flag(""); // remove -h,--help, this function is only to check priority!
    tmp.allow_extras();
    std::string str;
    addPriorityParam(tmp, str);
    if(argc-1)
    {
        try
        {
            (tmp).parse((argc-1), (argv+1));
        } 
        catch (const CLI::ParseError& e)
        {        
            std::cerr << "Internal error\n"; //should never happen
            int c = (tmp).exit(e);
            exit(c);
        }
    }
    return compressionPriorityFromString(str);
}

std::string commandLineToString(int argc, char** argv)
{
    std::string res;
    for (int i = 0; i < argc - 1; ++i)
        res += argv[i] + std::string(" ");
    res += argv[argc - 1];
    return res;
}

int parse_params(int argc, char** argv)
{
    CInfo info;
    info.time = get_time();
    info.full_command_line = commandLineToString(argc, argv);
    
	CLI::App app{ std::string("CoLoRd: Compressing long reads\n") +
        "version: " + std::to_string(version_major) + "." + std::to_string(version_minor) + "." + std::to_string(version_patch)
    };
    
    CDecompressorParams decomParams;

    CCompressorParams comParamsONT, comParamsPBRaw, comParamsPBHiFi;

    switch (checkPriority(argc, argv))
    {
    case CompressionPriority::Ratio:
        comParamsONT = compr_ONT_ratio_set_defaults();
        comParamsPBRaw = compr_PBRaw_ratio_set_defaults();        
        comParamsPBHiFi = compr_PBHiFi_ratio_set_defaults();
        break;
    case CompressionPriority::Balanced:
        comParamsONT = compr_ONT_balanced_set_defaults();
        comParamsPBRaw = compr_PBRaw_balanced_set_defaults();        
        comParamsPBHiFi = compr_PBHiFi_balanced_set_defaults();
        break;
    case CompressionPriority::Memory:
        comParamsONT = compr_ONT_memory_set_defaults();
        comParamsPBRaw = compr_PBRaw_memory_set_defaults();
        comParamsPBHiFi = compr_PBHiFi_memory_set_defaults();
        break;    
    }

    std::vector<CLI::Option*> toHideIfNoHelp;

    set_compression_params(&app, info, "compress-ont", "compress ONT data", comParamsONT, false, toHideIfNoHelp);
    set_compression_params(&app, info, "compress-pbraw", "compress PacBio Raw data", comParamsPBRaw, false, toHideIfNoHelp);
    set_compression_params(&app, info, "compress-pbhifi", "compress PacBio HiFi data", comParamsPBHiFi, false, toHideIfNoHelp);

    CLI::App* decompParser = app.add_subcommand("decompress", "decompression mode");

    decompParser->callback([&decomParams]() {
        runDecompression(decomParams);
    });

    //positionals    
    decompParser->add_option("input", decomParams.inputFilePath, "archive path")->required(true)->check(CLI::ExistingFile);
    decompParser->add_option("output", decomParams.outputFilePath, "output file path")->required(true);

    decompParser->add_option("-G,--reference-genome", decomParams.refGenomePath,
        "optional reference genome path (multi-FASTA gzipped or not), required for reference-based archives with no reference genome embedded (`-G` compression without `-s` switch)")->check(CLI::ExistingFile);

    decompParser->add_flag("-v,--verbose", decomParams.verbose, "verbose");

    CInfoParams infoParams;
    CLI::App* infoParser = app.add_subcommand("info", "print archive informations");
    infoParser->add_option("input", infoParams.inputFilePath, "archive path")->required(true)->check(CLI::ExistingFile);

    infoParser->callback([&infoParams] {
        runInfo(infoParams);
    });

    app.prefix_command(true);
    app.require_subcommand();

    try 
    {
        (app).parse((argc), (argv));
    }
    catch (const CLI::CallForHelp& e)
    {        
        return (app).exit(e);
    }
    catch (const CLI::ParseError& e)
    {
        for (auto ptr : toHideIfNoHelp)
            ptr->group("");
        std::cout << app.help();
        return (app).exit(e);
    }

    return 0;	
}