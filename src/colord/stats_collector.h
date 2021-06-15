#pragma once
#include <cinttypes>
#include <limits>
#include <vector>
#include <string_view>

class CStatsCollector
{
	friend class CGlobalStatsCollector;
	struct StatsDetail
	{
		struct CCompressionStats
		{
			uint64_t n_alternative_left_flank{};
			uint64_t n_alternative_right_flank{};
			uint64_t n_alternative_in_between{};

			uint64_t n_plain_symbols{};

			uint64_t n_symb_coded_with_edit_script{};

			uint64_t n_edit_script_symbols{};

			uint64_t n_substitution{};
			uint64_t n_match{};
			uint64_t n_insertion{};
			uint64_t n_deletion{};

			uint64_t n_symb_anchors{};
			uint64_t n_anchors{};

			uint64_t n_left_flank_symb{};
			uint64_t n_right_flank_symb{};

		};
		//reads
		uint64_t n_reads = 0;
		uint64_t min_read_len = std::numeric_limits<uint32_t>::max();
		uint64_t max_read_len = 0ul;
		uint64_t tot_read_len = 0;

		//refuse reasons
		uint32_t n_not_enough_unique_mmers_in_enc_read{};
		uint32_t n_too_many_matches{};
		uint32_t n_too_low_anchors{};

		//compression stats
		std::vector<CCompressionStats> compression_stats;
		uint64_t n_plain_reads_tot{}; //n_plain_tot
		uint64_t n_plain_symb{};

		uint64_t n_plain_reads_with_n_tot{};
		uint64_t n_plain_with_n_symb{};

		uint32_t n_non_rev_choosen{};
		uint32_t n_rev_choosen{};
	};

	StatsDetail stats;
	bool verbose;
public:

	CStatsCollector(bool verbose) :
		verbose(verbose) 
	{
	}

	void LogNonRevChoosen()
	{
		++stats.n_non_rev_choosen;
	}
	void LogRevChoosen()
	{
		++stats.n_rev_choosen;
	}

	void LogRead(uint64_t len);
	void LogPlainReadWithN(uint64_t len)
	{
		++stats.n_plain_reads_with_n_tot;
		stats.n_plain_with_n_symb += len;
	}
	void LogPlainRead(uint64_t len)
	{
		++stats.n_plain_reads_tot;
		stats.n_plain_symb += len;
	}
	void LogLevel(uint32_t level)
	{
		if (stats.compression_stats.size() < level + 1ul)
			stats.compression_stats.resize(level + 1ul);
	}
	StatsDetail::CCompressionStats& ComprStats(uint32_t level)
	{
		return stats.compression_stats[level];
	}
	void LogEditScript(uint32_t level, std::string_view editScript, uint32_t n_enc_read_symb);
	void LogRefuse_not_enough_unique_mmers_in_enc_read()
	{
		++stats.n_not_enough_unique_mmers_in_enc_read;
	}
	void LogRefuse_to_many_matches()
	{
		++stats.n_too_many_matches;
	}
	void LogRefuse_too_low_anchors()
	{
		++stats.n_too_low_anchors;
	}
	~CStatsCollector();
};
