#include "stats_collector.h"
#include "utils.h"
#include <iostream>
#include <mutex>

#define USE_LogEditScript

class CGlobalStatsCollector
{
	std::mutex mtx;
	CStatsCollector::StatsDetail stats;
	bool print = false;
public:
	static CGlobalStatsCollector& Inst()
	{
		static CGlobalStatsCollector inst;
		return inst;
	}
	void Log(CStatsCollector& thread_stats)
	{
		std::lock_guard lck(mtx);
		print = true;
		if (stats.compression_stats.size() < thread_stats.stats.compression_stats.size())
			stats.compression_stats.resize(thread_stats.stats.compression_stats.size());

		for (uint32_t i = 0; i < thread_stats.stats.compression_stats.size(); ++i)
		{
			auto& dest = stats.compression_stats[i];
			const auto& src = thread_stats.stats.compression_stats[i];
			

			dest.n_alternative_left_flank += src.n_alternative_left_flank;
			dest.n_alternative_right_flank += src.n_alternative_right_flank;
			dest.n_alternative_in_between += src.n_alternative_in_between;

			dest.n_plain_symbols += src.n_plain_symbols;

			dest.n_symb_coded_with_edit_script += src.n_symb_coded_with_edit_script;

			dest.n_edit_script_symbols += src.n_edit_script_symbols;

			dest.n_substitution += src.n_substitution;
			dest.n_match += src.n_match;
			dest.n_insertion += src.n_insertion;
			dest.n_deletion += src.n_deletion;

			dest.n_symb_anchors += src.n_symb_anchors;
			dest.n_anchors += src.n_anchors;

			dest.n_left_flank_symb += src.n_left_flank_symb;
			dest.n_right_flank_symb += src.n_right_flank_symb;
		}

		//reads
		stats.n_reads += thread_stats.stats.n_reads;
		if (thread_stats.stats.min_read_len < stats.min_read_len)
			stats.min_read_len = thread_stats.stats.min_read_len;
		if (thread_stats.stats.max_read_len > stats.max_read_len)
			stats.max_read_len = thread_stats.stats.max_read_len;

		stats.tot_read_len += thread_stats.stats.tot_read_len;

		//refuse reasons
		stats.n_not_enough_unique_mmers_in_enc_read += thread_stats.stats.n_not_enough_unique_mmers_in_enc_read;
		stats.n_too_many_matches += thread_stats.stats.n_too_many_matches;
		stats.n_too_low_anchors += thread_stats.stats.n_too_low_anchors;

		//compression stats
		
		stats.n_plain_reads_tot += thread_stats.stats.n_plain_reads_tot;
		stats.n_plain_symb += thread_stats.stats.n_plain_symb;

		stats.n_plain_reads_with_n_tot += thread_stats.stats.n_plain_reads_with_n_tot;
		stats.n_plain_with_n_symb += thread_stats.stats.n_plain_with_n_symb;

		stats.n_non_rev_choosen += thread_stats.stats.n_non_rev_choosen;
		stats.n_rev_choosen += thread_stats.stats.n_rev_choosen;
	}

	~CGlobalStatsCollector()
	{
		if (!print)
			return;
		std::ostream& summary = std::cerr;

		summary << " * * * * * * * * READS STATS * * * * * * * * \n";
		summary << "# reads                        : " << stats.n_reads << "\n";
		summary << "min read len                   : " << stats.min_read_len << "\n";
		summary << "max read len                   : " << stats.max_read_len << "\n";
		summary << "# symbols                      : " << stats.tot_read_len << "\n";

		summary << " * * * * * * * * REFUSE REASONS STATS * * * * * * * * \n";
		summary << "# not enough uniq mmers in enc : " << stats.n_not_enough_unique_mmers_in_enc_read << "\n";
		summary << "# too many matches             : " << stats.n_too_many_matches << "\n";
		summary << "# too low anchors              : " << stats.n_too_low_anchors << "\n";
		
		summary << " * * * * * * * * COMPRESSION STATS * * * * * * * * \n";
		summary << "# plain reads                  : " << stats.n_plain_reads_tot << "\n";
		summary << "# symb plain reads             : " << stats.n_plain_symb << "\n";

		summary << "# plain reads (reason: N)      : " << stats.n_plain_reads_with_n_tot << "\n";
		summary << "# symb plain reads (reason: N) : " << stats.n_plain_with_n_symb << "\n";

		summary << "# non rev choosen              : " << stats.n_non_rev_choosen << "\n";
		summary << "# rev choosen                  : " << stats.n_rev_choosen << "\n";

		for (uint32_t i = 0; i < stats.compression_stats.size(); ++i)
		{
			summary << " --------------- level " << i << " --------------- \n";		
			summary << "# alt for left flank           : " << stats.compression_stats[i].n_alternative_left_flank << "\n";
			summary << "# alt in between anchors       : " << stats.compression_stats[i].n_alternative_in_between << "\n";
			summary << "# alt for right flank          : " << stats.compression_stats[i].n_alternative_right_flank << "\n";					
			summary << "# symb plain                   : " << stats.compression_stats[i].n_plain_symbols << "\n";
			summary << "# symb edit script encoded     : " << stats.compression_stats[i].n_symb_coded_with_edit_script << "\n";
			summary << "# symb in edit script          : " << stats.compression_stats[i].n_edit_script_symbols << "\n";
			summary << "# mismatches                   : " << stats.compression_stats[i].n_substitution << "\n";
			summary << "# matches                      : " << stats.compression_stats[i].n_match << "\n";
			summary << "# insertions                   : " << stats.compression_stats[i].n_insertion << "\n";
			summary << "# deletions                    : " << stats.compression_stats[i].n_deletion << "\n";
			summary << "# symb anchors                 : " << stats.compression_stats[i].n_symb_anchors << "\n";
			summary << "# anchors                      : " << stats.compression_stats[i].n_anchors << "\n";
			summary << "# symb left flank              : " << stats.compression_stats[i].n_left_flank_symb << "\n";
			summary << "# symb right flank             : " << stats.compression_stats[i].n_right_flank_symb << "\n";
		}
	}
};

CStatsCollector::~CStatsCollector()
{
	if(verbose)
		CGlobalStatsCollector::Inst().Log(*this);
}

void CStatsCollector::LogRead(uint64_t len)
{
	stats.tot_read_len += len;
	if (len > stats.max_read_len)
		stats.max_read_len = len;
	else if (len < stats.min_read_len)
		stats.min_read_len = len;
	++stats.n_reads;
}

void CStatsCollector::LogEditScript(uint32_t level, std::string_view editScript, uint32_t n_enc_read_symb)
{
#ifdef USE_LogEditScript
	stats.compression_stats[level].n_symb_coded_with_edit_script += n_enc_read_symb;
	stats.compression_stats[level].n_edit_script_symbols += editScript.length();
	uint64_t n_subs{}, n_ins{}, n_del{}, n_match{};
	for(auto s : editScript)
		if (CMissmatchCoder::is_missmatch_es_symb(s))
			++n_subs;
		else if (s == 'A' || s == 'C' || s == 'G' || s == 'T')
			++n_ins;
		else if (s == 'D')
			++n_del;
		else if (s == 'M')
			++n_match;
		else
		{
			std::cerr << "Error: wrong symb in edit script: " << s << "\n";
			exit(1);
		}
	stats.compression_stats[level].n_substitution += n_subs;
	stats.compression_stats[level].n_match += n_match;
	stats.compression_stats[level].n_insertion += n_ins;
	stats.compression_stats[level].n_deletion += n_del;
#endif
}