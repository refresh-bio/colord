#pragma once

#include <cctype>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

inline uint64_t get_time()
{
	return static_cast<uint64_t>(time(nullptr));
}


enum class QualityComprMode {
	Original,
	QuinaryAverage,
	QuadAverage,
	BinaryAverage,
	QuinaryThreshold,
	QuadThreshold,
	BinaryThreshold,
	Average,
	None
};
enum class HeaderComprMode { Original, Main, None };
enum class ReferenceReadsMode { All, Sparse }; //All - all reads (except containing N) are kept as reference, Sparse - only some part of reads is kept as reference
enum class DataSource {ONT, PBRaw, PBHiFi};
enum class CompressionPriority { Ratio, Balanced, Memory };
struct CCompressorParams
{
	DataSource dataSource = DataSource::ONT;
	CompressionPriority priority = CompressionPriority::Balanced;

	std::string inputFilePath;
	std::string outputFilePath;
	uint32_t kmerLen = 0;
	uint32_t anchorLen = 0;
	int32_t compressionLevel = 2;
	uint32_t minKmerCount = 3;
	uint32_t maxKmerCount = 100;
	uint32_t filterHashModulo = 9;
	uint32_t maxCandidates = 8;
	double editScriptCostMultiplier = 1.0;
	uint32_t maxRecurence = 5;
	uint32_t minPartLenToConsiderAltRead = 48; // minimum length of encoding part to consider using alternative read
	double minFractionOfMmersInEncode = 0.5; // if A is set of m-mers in encode read R then read is refused from encoding if |A| < minFractionOfMmersInEncode * len(R)
	double minFractionOfMmersInEncodeToAlwaysEncode = 0.9; // if A is set of m-mers in encode read R then read is accepted to encoding always if |A| > minFractionOfMmersInEncodeToAlwaysEncode * len(R)
	double maxMatchesMultiplier = 10; // if the number of matches between encode read R and reference read is r, then read is refused from encodinf if r > maxMatchesMultiplier * len(R)
	uint32_t nThreads = std::max(std::thread::hardware_concurrency(), 1u);
	QualityComprMode qualityComprMode = QualityComprMode::QuadThreshold;
	//int qualityThreshold = 7;
	std::vector<uint32_t> qualityFwdThresholds;
	std::vector<uint32_t> qualityRevThresholds;
	HeaderComprMode headerComprMode = HeaderComprMode::Original;
	ReferenceReadsMode referenceReadsMode = ReferenceReadsMode::Sparse;
	//uint32_t sparseMode_range_symbols = 400; // only for sparse reference reads mode, in million of symbols
	double sparseMode_range_symbols = 2; // only for sparse reference reads mode, the number of reference reads in range is determined based on the number of symbols, which in turn is determined by the number of trusted unique k-mers multiplied by the value of this parameter
	double sparseMode_exponent = 1.0; // only for sparse reference reads mode

	uint32_t minAnchors = 1; //if number of anchors common to encode read and reference candidate is lower than minAnchors candidate is refused

	bool verbose = false;

	double fillFactorFilteredKmers = 0.75;
	double fillFactorKmersToReads = 0.8;

	std::string refGenomePath;
	bool storeRefGenome = false;

	struct {
		std::string qual_mode;
		std::string header_mode;
		std::string reference_mode;
		std::string priority;
	} internal; //for parsing
};

CompressionPriority compressionPriorityFromString(const std::string& str);
std::string compressionPriorityToString(CompressionPriority compressionPriority);

std::string dataSourceToString(DataSource dataSource);

HeaderComprMode headerComprModeFromString(const std::string& str);
std::string headerComprModeToString(HeaderComprMode headerComprMode);

std::string qualityComprModeToString(QualityComprMode qualityComprMode);

ReferenceReadsMode referenceReadsModeFromString(const std::string& str);
std::string referenceReadsModeToString(ReferenceReadsMode referenceReadsMode);

std::string vec_to_string(const std::vector<uint32_t>& vec);


struct CDecompressorParams
{
	std::string inputFilePath;
	std::string outputFilePath;

	std::string refGenomePath;
	bool verbose = false;
};

struct CInfoParams
{
	std::string inputFilePath;
};