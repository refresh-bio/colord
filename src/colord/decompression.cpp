#include "decompression.h"
#include "utils.h"
#include "reference_genome.h"
#include "dna_coder.h"
#include "archive.h"
#include "parallel_queue.h"
#include "entr_read.h"
#include "entr_header.h"
#include "entr_qual.h"
#include <cstring>

// *****************************************************************
void store_read(ofstream& out, read_t& read)
{
	for (uint32_t i = 0; i < read.size() - 1; ++i)
		read[i] = "ACGTN"[read[i]];

	read.back() = '\n';

	out.write((char*)read.data(), read.size());
}

class CBufferedWrite
{
	std::ofstream out;
	const size_t write_buff_size = 1ull << 26;
	std::vector<uint8_t> write_buf;
	size_t pos{};
public:
	CBufferedWrite(const std::string& outputFilePath) :
		out(outOpenOrDie(outputFilePath, std::ios::binary)),
		write_buf(write_buff_size)
	{

	}
	void write(const char* data, size_t size)
	{		
		while (size)
		{
			if (pos + size <= write_buff_size)
			{
				memcpy(write_buf.data() + pos, data, size);
				pos += size;
				return;
			}
			else
			{
				size_t part_size = write_buff_size - pos;
				memcpy(write_buf.data() + pos, data, part_size);
				pos += part_size;
				out.write(reinterpret_cast<char*>(write_buf.data()), pos);
				data += part_size;				
				size -= part_size;				
				pos = 0;
			}
		}
	}
	~CBufferedWrite()
	{
		out.write(reinterpret_cast<char*>(write_buf.data()), pos);
	}
};

class CFastaWriter
{
	CBufferedWrite out;
	CParallelQueue<decomp_read_pack_t>& read_decompr_queue;
	CParallelQueue<header_pack_t>& header_decompr_queue;	

	decomp_read_pack_t read_pack;
	uint32_t read_pack_pos{};

	header_pack_t header_pack;
	uint32_t header_pack_pos{};

	
	template<typename T_ELEM, typename PACK_T>
	bool nextImpl(T_ELEM& elem, CParallelQueue<PACK_T>& queue, PACK_T& pack, uint32_t& pos)
	{
		if (pos == pack.size())
		{
			pos = 0;
			if (!queue.Pop(pack))
				return false;
		}
		elem = std::move(pack[pos++]);
		return true;
	}

	bool nextRead(read_t& read)
	{
		return nextImpl(read, read_decompr_queue, read_pack, read_pack_pos);
	}
	
	bool nextHeader(header_elem_t& header)
	{
		return nextImpl(header, header_decompr_queue, header_pack, header_pack_pos);
	}
public:
	CFastaWriter(const std::string& outputFilePath,
		CParallelQueue<decomp_read_pack_t>& read_decompr_queue,
		CParallelQueue<header_pack_t>& header_decompr_queue
		) :
		out(outputFilePath),
		read_decompr_queue(read_decompr_queue),
		header_decompr_queue(header_decompr_queue)		
	{

	}

	void Write()
	{
		read_t read;
		qual_t qual;
		header_elem_t header;

		bool hasNextRead = nextRead(read);		
		bool hasNextHeader = nextHeader(header);
		while (hasNextRead && hasNextHeader)
		{
			out.write(">", 1);
			out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);

			for (uint32_t i = 0; i < read.size() - 1; ++i)
				read[i] = "ACGTN"[read[i]];
			read.back() = '\n';
			out.write((char*)read.data(), read.size());
			
			hasNextRead = nextRead(read);			
			hasNextHeader = nextHeader(header);
		}
		if (hasNextRead + hasNextHeader != 0)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}
	}
};

class CFastqWriter
{
	CBufferedWrite out;
	CParallelQueue<decomp_read_pack_t>& read_decompr_queue;
	CParallelQueue<header_pack_t>& header_decompr_queue;
	CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue;

	
	decomp_read_pack_t read_pack;
	uint32_t read_pack_pos{};

	header_pack_t header_pack;
	uint32_t header_pack_pos{};

	decomp_qual_pack_t qual_pack;
	uint32_t qual_pack_pos{};
	

	template<typename T_ELEM, typename PACK_T>
	bool nextImpl(T_ELEM& elem, CParallelQueue<PACK_T>& queue, PACK_T& pack, uint32_t& pos)
	{
		if (pos == pack.size())
		{
			pos = 0;
			if (!queue.Pop(pack))
				return false;
		}
		elem = std::move(pack[pos++]);
		return true;
	}

	bool nextRead(read_t& read)
	{
		return nextImpl(read, read_decompr_queue, read_pack, read_pack_pos);		
	}
	bool nextQual(qual_t& qual)
	{
		return nextImpl(qual, qual_decompr_queue, qual_pack, qual_pack_pos);
	}
	bool nextHeader(header_elem_t& header)
	{
		return nextImpl(header, header_decompr_queue, header_pack, header_pack_pos);
	}
	void write_via_buff(uint8_t* data, size_t size)
	{

	}
public:
	CFastqWriter(const std::string& outputFilePath,
		CParallelQueue<decomp_read_pack_t>& read_decompr_queue,
		CParallelQueue<header_pack_t>& header_decompr_queue,
		CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue) :
		out(outputFilePath),
		read_decompr_queue(read_decompr_queue),
		header_decompr_queue(header_decompr_queue),
		qual_decompr_queue(qual_decompr_queue)
		
	{

	}

	void Write()
	{
		read_t read;
		qual_t qual;
		header_elem_t header;

		bool hasNextRead = nextRead(read);
		bool hasNextQual = nextQual(qual);
		bool hasNextHeader = nextHeader(header);
		while (hasNextRead && hasNextQual && hasNextHeader)
		{
			out.write("@", 1);
			out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);

			for (uint32_t i = 0; i < read.size()-1; ++i)
				read[i] = read[i] = "ACGTN"[read[i]];
			read.back() = '\n';
			out.write((char*)read.data(), read.size());

			out.write("+", 1);
			if(header.second == qual_header_type::eq_read_header)
				out.write(header.first.c_str(), header.first.size());
			out.write("\n", 1);
			
//			out.write(qual.c_str(), qual.size());
			out.write((char*)qual.data(), qual.size());
			out.write("\n", 1);

			hasNextRead = nextRead(read);
			hasNextQual = nextQual(qual);
			hasNextHeader = nextHeader(header);
		}
		if (hasNextRead + hasNextQual + hasNextHeader != 0)
		{
			std::cerr << "Error: cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__ << "\n";
			exit(1);
		}
	}
};

bool checkVersion(CInfo& info)
{
	return info.version_major == version_major;	
}

void runDecompression(const CDecompressorParams& params)
{
	std::cerr << "Running decompression.\n";
	CArchive archive(true);
	if (!archive.Open(params.inputFilePath))
	{
		std::cerr << "Error: cannot open archive: " << params.inputFilePath << "\n";
		exit(1);
	}
	CInfo info;
	int s_info = archive.GetStreamId("info");
	size_t meta;
	std::vector<uint8_t> _info;
	archive.GetPart(s_info, _info, meta);
	info.Deserialize(_info);


	if (!checkVersion(info))
	{
		std::cerr << "Error: incompatibile archive version\n";
		exit(1);
	}

	int s_qual = archive.GetStreamId("qual");
	int s_meta = archive.GetStreamId("meta");

	bool is_fastq = s_qual != -1;

	uint32_t maxCandidates{};
	int32_t compressionLevel{};
	DataSource dataSource{};
	uint64_t approx_stream_size; //approximate size of DNA (or qual, cause the are the same) BEFORE compression (or after decompression)

	std::vector<uint8_t> _params;
	
	uint32_t tot_ref_reads{};

	archive.GetPart(s_meta, _params, meta);
	const uint8_t* ptr_params = _params.data();

	LoadLittleEndian(ptr_params, tot_ref_reads);
	ptr_params += sizeof(tot_ref_reads);

	LoadLittleEndian(ptr_params, maxCandidates);
	ptr_params += sizeof(maxCandidates);

	LoadLittleEndian(ptr_params, compressionLevel);
	ptr_params += sizeof(compressionLevel);

	uint8_t tmp;
	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	dataSource = static_cast<DataSource>(tmp);

	LoadLittleEndian(ptr_params, approx_stream_size);
	ptr_params += sizeof(approx_stream_size);

	QualityComprMode qualityComprMode{};
	vector<uint32_t> qualityFwdThresholds;
	vector<uint32_t> qualityRevThresholds;
	if (is_fastq)
	{
		LoadLittleEndian(ptr_params, tmp);
		ptr_params += sizeof(tmp);
		qualityComprMode = static_cast<QualityComprMode>(tmp);

		switch (qualityComprMode)
		{
		case QualityComprMode::None:
			qualityRevThresholds.emplace_back();
			LoadLittleEndian(ptr_params, qualityRevThresholds.back());
			ptr_params += sizeof(qualityRevThresholds.back());
			break;
		case QualityComprMode::Original:
			break;
		case QualityComprMode::BinaryThreshold:
			for (uint32_t i = 0; i < 2; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::QuadThreshold:
			for (uint32_t i = 0; i < 4; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::QuinaryThreshold:
			for (uint32_t i = 0; i < 5; ++i)
			{
				qualityRevThresholds.emplace_back();
				LoadLittleEndian(ptr_params, qualityRevThresholds.back());
				ptr_params += sizeof(qualityRevThresholds.back());
			}
			break;
		case QualityComprMode::Average:
			break;
		case QualityComprMode::BinaryAverage:
			break;
		case QualityComprMode::QuadAverage:
			break;
		case QualityComprMode::QuinaryAverage:
			break;
		default:
			break;
		}

		if(params.verbose)
		{
			std::cerr << "qualityRevThresholds: ";
			for (auto v : qualityRevThresholds)
				std::cerr << v << " ";
			std::cerr << "\n";
			switch (qualityComprMode)
			{
			case QualityComprMode::None:
				std::cerr << "qualityComprMode: None\n";
				break;

			case QualityComprMode::Original:
				std::cerr << "qualityComprMode: Original\n";
				break;
			case QualityComprMode::BinaryThreshold:
				std::cerr << "qualityComprMode: BinaryThreshold\n";
				break;
			case QualityComprMode::QuadThreshold:
				std::cerr << "qualityComprMode: QuadThreshold\n";
				break;
			case QualityComprMode::QuinaryThreshold:
				std::cerr << "qualityComprMode: QuinaryThreshold\n";
				break;
			case QualityComprMode::Average:
				std::cerr << "qualityComprMode: Average\n";
				break;
			case QualityComprMode::BinaryAverage:
				std::cerr << "qualityComprMode: BinaryAverage\n";
				break;
			case QualityComprMode::QuadAverage:
				std::cerr << "qualityComprMode: QuadAverage\n";
				break;
			case QualityComprMode::QuinaryAverage:
				std::cerr << "qualityComprMode: QuinaryAverage\n";
				break;
			default:
				break;
			}
		}
	}
	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	HeaderComprMode headerComprMode = static_cast<HeaderComprMode>(tmp);

	if (params.verbose)
	{
		switch (headerComprMode)
		{
		case HeaderComprMode::Original:
			std::cerr << "headerComprMode: Original\n";
			break;
		case HeaderComprMode::Main:
			std::cerr << "headerComprMode: Main\n";
			break;
		case HeaderComprMode::None:
			std::cerr << "headerComprMode: None\n";
			break;
		default:
			break;
		}
	}

	LoadLittleEndian(ptr_params, tmp);
	ptr_params += sizeof(tmp);
	ReferenceReadsMode referenceReadsMode = static_cast<ReferenceReadsMode>(tmp);

	if (params.verbose)
	{
		switch (referenceReadsMode)
		{
		case ReferenceReadsMode::All:
			std::cerr << "reference reads mode: All\n";
			break;
		case ReferenceReadsMode::Sparse:
			std::cerr << "reference reads mode: Sparse\n";
			break;
		default:
			break;
		}
	}

	uint32_t sparseRange = 0;
	double sparseExponent = 0;
	if (referenceReadsMode == ReferenceReadsMode::Sparse)
	{
		LoadLittleEndian(ptr_params, sparseRange);
		ptr_params += sizeof(sparseRange);
		
		LoadLittleEndian(ptr_params, sparseExponent);
		ptr_params += sizeof(sparseExponent);
		if(params.verbose)
		{
			std::cerr << "sparse range: " << sparseRange << "\n";
			std::cerr << "sparse exponent: " << sparseExponent << "\n";
		}
	}

	bool ref_genome_available;
	LoadLittleEndian(ptr_params, ref_genome_available);
	ptr_params += sizeof(ref_genome_available);

	bool ref_genome_in_arch;
	
	uint32_t n_ref_genome_pseudo_reads{};
	
	std::unique_ptr<CReferenceGenome> ref_genome;
	if (ref_genome_available)
	{
		std::vector<uint8_t> ref_genome_checksum;
		//uint64_t mean_read_len{};
		//uint32_t kmerLen{};

		uint32_t ref_genome_overlap_size;
		uint32_t ref_genome_read_len;

		LoadLittleEndian(ptr_params, ref_genome_in_arch);
		ptr_params += sizeof(ref_genome_in_arch);

		LoadLittleEndian(ptr_params, ref_genome_read_len);
		ptr_params += sizeof(ref_genome_read_len);

		LoadLittleEndian(ptr_params, ref_genome_overlap_size);
		ptr_params += sizeof(ref_genome_overlap_size);

		LoadLittleEndian(ptr_params, n_ref_genome_pseudo_reads);
		ptr_params += sizeof(n_ref_genome_pseudo_reads);

		if (!ref_genome_in_arch)
		{
			ref_genome_checksum.resize(16);
			for(auto& c: ref_genome_checksum)
			{
				LoadLittleEndian(ptr_params, c);
				ptr_params += sizeof(c);			
			}
		}

		
		if (ref_genome_in_arch)
			ref_genome = std::make_unique<CReferenceGenome>(archive, ref_genome_overlap_size, params.verbose);
		else
		{
			if (params.refGenomePath == "")
			{
				std::cerr << "Error: compressed file was created without -s switch, reference genome is required for decompression";
				exit(1);
			}
			ref_genome = std::make_unique<CReferenceGenome>(params.refGenomePath, ref_genome_overlap_size, true, params.verbose);
			if (ref_genome->GetChecksum() != ref_genome_checksum)
			{
				std::cerr << "Error: different reference genome was used during compression. Decompression impossible.\n";
				exit(1);
			}
		}

		ref_genome->SetReadLen(ref_genome_read_len);
	}

	CReferenceReads ref_reads(tot_ref_reads);

	if (ref_genome_available)
	{
		read_pack_t reads = ref_genome->GetPseudoReads();
		for (const auto& [tmp, read] : reads)
			ref_reads.Add(read);
	}

	CRefReadsAccepter ref_reads_accepter(sparseRange, sparseExponent, n_ref_genome_pseudo_reads);

	if(params.verbose)
	{
		std::cerr << "maxCandidates: " << maxCandidates << "\n";
		std::cerr << "compressionLevel: " << compressionLevel << "\n";
		std::cerr << "approx stream size: " << approx_stream_size << "\n";
	}

	std::vector<CParallelQueue<decomp_read_pack_t>*> read_decompr_queues;

	CParallelQueue<decomp_read_pack_t> read_decompr_queue(read_decompress_queue_size);
	read_decompr_queues.push_back(&read_decompr_queue);

	std::unique_ptr<CParallelQueue<decomp_read_pack_t>> read_decompr_queue_for_qual;
	std::unique_ptr<CParallelQueue<decomp_qual_pack_t>> qual_decompr_queue;
	CParallelQueue<header_pack_t> header_decompr_queue(header_decompress_queue_size);

	if (is_fastq)//second queue for quality decompression
	{
		read_decompr_queue_for_qual = std::make_unique<CParallelQueue<decomp_read_pack_t>>(read_decompress_queue_size);
		read_decompr_queues.push_back(read_decompr_queue_for_qual.get());
		qual_decompr_queue = std::make_unique<CParallelQueue<decomp_qual_pack_t>>(qual_decompress_queue_size);
	}

	std::thread reads_decompr([&archive, &read_decompr_queues, &params, maxCandidates, referenceReadsMode, &ref_reads, &ref_reads_accepter, compressionLevel, approx_stream_size, n_ref_genome_pseudo_reads, &info] {
		CEntropyDecomprReads decompr{ archive, read_decompr_queues, params.verbose, maxCandidates, compressionLevel, approx_stream_size, info.total_reads, referenceReadsMode, ref_reads, ref_reads_accepter, n_ref_genome_pseudo_reads };
		decompr.Decompress();
	});
	
	
	std::thread headers_decompr([&archive, &header_decompr_queue, headerComprMode, compressionLevel, &params] {
		CEntrDecomprHeaders decompr{ archive, header_decompr_queue, headerComprMode, compressionLevel, params.verbose };
		decompr.Decompress();
	});
	
	std::thread qual_decompr;
	if (is_fastq)
	{
		qual_decompr = std::thread([&archive, &read_decompr_queues, &qual_decompr_queue, &params, qualityComprMode, dataSource, &qualityFwdThresholds, &qualityRevThresholds, compressionLevel, approx_stream_size] {
			CEntrDecomprQuals decompr{ archive, *read_decompr_queues[1], *qual_decompr_queue.get(), params.verbose, qualityComprMode, dataSource, qualityFwdThresholds, qualityRevThresholds, compressionLevel, approx_stream_size };
			decompr.Decompress();
		});
	}

	auto out = outOpenOrDie(params.outputFilePath, std::ios::binary);
		
	if (is_fastq)
	{
		CFastqWriter{ params.outputFilePath, read_decompr_queue, header_decompr_queue, *qual_decompr_queue.get() }.Write();
		qual_decompr.join();
	}
	else
	{
		CFastaWriter{ params.outputFilePath, read_decompr_queue, header_decompr_queue }.Write();
	}

	headers_decompr.join();
	reads_decompr.join();

	archive.Close();

	//std::cerr << "n reads: " << n_reads << "\n";
	//coder.Finish();
}
