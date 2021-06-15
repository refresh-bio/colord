#include "reads_sim_graph.h"
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <map>
#include <iterator>

#include "pooled_threads.h"
#include "timer.h"


#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
class GetAcceptedKmersQueue
{
	struct Elem
	{
		uint32_t start, end;
	};
	std::queue<Elem> q;
	std::mutex mtx;
	std::condition_variable cv_pop;
	std::condition_variable cv_sync;
	uint32_t allowed, running{};
	bool finished = false;
public:
	GetAcceptedKmersQueue()
	{

	}
	void Push(uint32_t start, uint32_t end)
	{
		q.emplace(Elem{ start, end });
	}

	//call at the beginning
	void SetAllowed(uint32_t allowed)
	{
		std::lock_guard<std::mutex> lck(mtx);
		this->allowed = allowed;
	}

	bool Pop(uint32_t& start, uint32_t& end)
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.wait(lck, [this] {
			return finished || (!q.empty() && running < allowed);
		});
		if (finished && q.empty())
			return false;
		start = q.front().start;
		end = q.front().end;
		q.pop();
		++running;
		return true;
	}

	void NotifyTaskCompleted(uint32_t newAllowed)
	{
		std::lock_guard<std::mutex> lck(mtx);
		--running;
		allowed = newAllowed;
		if (q.empty() && !running) //it was last running threads, warning: it will work only if all task are produced before tasks starts
			cv_sync.notify_one();
		else
			cv_pop.notify_all();
	}
	void Finish()
	{
		std::lock_guard<std::mutex> lck(mtx);
		finished = true;
		cv_pop.notify_all();
	}

	void Run()
	{
		std::unique_lock<std::mutex> lck(mtx);
		cv_pop.notify_all();
		cv_sync.wait(lck);
	}
};

class CReadsSimilarityGraphInternalThreads
{
	std::vector<std::thread> threads;
#ifdef MEASURE_THREADS_TIMES
	std::vector<double> twTimes;
#endif
	GetAcceptedKmersQueue Q;
	uint32_t kmer_len;
	const CKmerFilter& filteredKmers;
	CParallelQueuePopWaiting<CCompressPack>& compress_queue;
	const read_pack_t* reads_pack;
	std::vector<std::pair<read_t*, std::vector<kmer_type>>>* accepted_kmers;
public:
#ifdef MEASURE_THREADS_TIMES
	std::vector<double>& GetTwTimes()
	{
		return twTimes;
	}
#endif
	CReadsSimilarityGraphInternalThreads(uint32_t maxRunning, uint32_t _kmer_len, const CKmerFilter& _filteredKmers, CParallelQueuePopWaiting<CCompressPack>& _compress_queue) :
		kmer_len(_kmer_len),
		filteredKmers(_filteredKmers),
		compress_queue(_compress_queue)
	{
#ifdef MEASURE_THREADS_TIMES
		twTimes.resize(maxRunning);
#endif
		for (uint32_t tno = 0; tno < maxRunning; ++tno)
			threads.emplace_back([this, tno] {
#ifdef MEASURE_THREADS_TIMES
			CThreadWatch tw;
			tw.startTimer();
#endif
			uint32_t s, e;
			while (Q.Pop(s, e))
			{
				for (uint32_t i = s; i < e; ++i)
				{
					auto& [hasN, read] = (*reads_pack)[i];
					(*accepted_kmers)[i].first = (read_t*)&read;

					if (hasN)
						continue;

					if (read_len(read) < kmer_len)
						continue;

					kmer_type kmer;
					CKmerWalker kmerWalker(read, kmer_len, kmer);
					uint32_t kmers_in_read = static_cast<uint32_t>(read_len(read)) - kmer_len + 1;
					hash_set_t already_processed_kmers(~0ull, static_cast<size_t>(kmers_in_read / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});

					while (kmerWalker.NextKmer())
					{
//						if (auto it = already_processed_kmers.insert(kmer); !it.second) //skip already processed k-mer
						if (!filteredKmers.Possible(kmer) || !already_processed_kmers.insert_fast(kmer)) //skip already processed k-mer
							continue;

						if (filteredKmers.Check(kmer))
							(*accepted_kmers)[i].second.push_back(kmer);
					}
				}
				Q.NotifyTaskCompleted(compress_queue.GetNWaitingOnPop() + 1);
			}
#ifdef MEASURE_THREADS_TIMES
			tw.stopTimer();
			twTimes[tno] += tw.getElapsedTime();
#endif
		});
	}

	void Run(const read_pack_t* reads_pack, std::vector<std::pair<read_t*, std::vector<kmer_type>>>* accepted_kmers)
	{
		this->reads_pack = reads_pack;
		this->accepted_kmers = accepted_kmers;
		Q.SetAllowed(compress_queue.GetNWaitingOnPop() + 1);

//		uint32_t approx_task_size = reads_pack->size() / (threads.size() + 1) / 8;
/*		uint32_t approx_task_size = reads_pack->size() / (threads.size() + 1) / 4;
		if (!approx_task_size)
			approx_task_size = 1;
		for (uint32_t start = 0; start < reads_pack->size(); )
		{
			uint32_t end = start + approx_task_size;
			if (end > reads_pack->size())
				end = reads_pack->size();
			Q.Push(start, end);
			start = end;
		}*/

		uint64_t n_parts = 3 * (threads.size() + 1);
		uint64_t init_task_size = 3 * reads_pack->size() / (2 * n_parts);
		if (!init_task_size)
			init_task_size = 1;
		double decreaser = (double) reads_pack->size() / n_parts / n_parts;
		double cur_task_size = static_cast<double>(init_task_size);

		for (uint64_t start = 0; start < reads_pack->size(); )
		{
			uint64_t end = start + (uint64_t) cur_task_size;
			if (end > static_cast<uint32_t>(reads_pack->size()))
				end = static_cast<uint32_t>(reads_pack->size());
			Q.Push(static_cast<uint32_t>(start), static_cast<uint32_t>(end));
			start = end;
			cur_task_size -= decreaser;
			if (cur_task_size < 1)
				cur_task_size = 1;
		}

		Q.Run();
	}
	void NotifyFinish()
	{
		Q.Finish();
		for (auto& t : threads)
			t.join();
	}
};


std::vector<std::pair<read_t*, std::vector<kmer_type>>> CReadsSimilarityGraph::getAcceptedKmers(const read_pack_t& reads_pack)
{
	std::vector<std::pair<read_t*, std::vector<kmer_type>>> accepted_kmers(reads_pack.size());

	internalThreads->Run(&reads_pack, &accepted_kmers);

	return accepted_kmers;
}
#else

std::vector<std::pair<read_t*, std::vector<kmer_type>>> CReadsSimilarityGraph::getAcceptedKmers(const read_pack_t& reads_pack)
{
	std::vector<std::pair<read_t*, std::vector<kmer_type>>> accepted_kmers(reads_pack.size());
	uint32_t n_threads = 1 + compress_queue.GetNWaitingOnPop();
	std::vector<pooled_threads::thread> threads;
	uint32_t per_threads = static_cast<uint32_t>(reads_pack.size()) / n_threads;

#ifdef MEASURE_THREADS_TIMES
	if(n_threads > twTimes.size())
		twTimes.resize(n_threads);
#endif

	for (uint32_t tno = 0; tno < n_threads; ++tno)
	{
		uint32_t s = tno * per_threads;
		uint32_t e = (tno + 1) * per_threads;
		if (tno == n_threads - 1)
			e = static_cast<uint32_t>(reads_pack.size());
		threads.emplace_back([&accepted_kmers, &reads_pack, this, s, e, tno] {
#ifdef MEASURE_THREADS_TIMES
			CThreadWatch tw;
			tw.startTimer();
#endif

			for (uint32_t i = s; i < e; ++i)				
			{
				auto& [hasN, read] = reads_pack[i];
				accepted_kmers[i].first = (read_t*)&read;

				if (hasN)
					continue;

				if (read_len(read) < kmer_len)
					continue;

				kmer_type kmer;
				CKmerWalker kmerWalker(read, kmer_len, kmer);
				uint32_t kmers_in_read = static_cast<uint32_t>(read_len(read)) - kmer_len + 1;
				hash_set_t already_processed_kmers(~0ull, static_cast<size_t>(kmers_in_read / 0.4), 0.4, std::equal_to<uint64_t>{}, MurMur64Hash{});

				while (kmerWalker.NextKmer())
				{
//					if (auto it = already_processed_kmers.insert(kmer); !it.second) //skip already processed k-mer
					if (!already_processed_kmers.insert_fast(kmer)) //skip already processed k-mer
						continue;

					if (filteredKmers.Check(kmer))
						accepted_kmers[i].second.push_back(kmer);
				}
			}
#ifdef MEASURE_THREADS_TIMES
			tw.stopTimer();
			twTimes[tno] += tw.getElapsedTime();
#endif
		});
	}

	for (auto& th : threads)
		th.join();

	return accepted_kmers;
}
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH

void CReadsSimilarityGraph::processReferenceGenome(CReferenceGenome* reference_genome)
{
	if (!reference_genome)
		return;

	read_pack_t reads_pack = reference_genome->GetPseudoReads();

	auto accepted_kmers = getAcceptedKmers(reads_pack);

	assert(reads_pack.size() == accepted_kmers.size());

	for (size_t i = 0; i < reads_pack.size(); ++i)
	{
		auto& [read, read_kmers] = accepted_kmers[i];
		
		reference_reads.Add(*read);
		++id_in_reference;
		
		uint32_t read_kmers_size = static_cast<uint32_t>(read_kmers.size());
		for (uint32_t i = 0; i < read_kmers_size; ++i)
		{			
			kmers.insert(read_kmers[i], id_in_reference - 1);
		}

		++current_read_id;
	}
	reference_genome->Release();
}

void CReadsSimilarityGraph::processReadsPack(const read_pack_t& reads_pack)
{
	auto accepted_kmers = getAcceptedKmers(reads_pack);

	std::vector<std::pair<uint32_t, uint32_t>> neighbours;
	std::unordered_map<uint32_t, uint32_t> m;

	m.max_load_factor(1.0);

	assert(reads_pack.size() == accepted_kmers.size());
	for(size_t i = 0 ; i < reads_pack.size() ; ++i)
	{
		auto& [read, read_kmers] = accepted_kmers[i];
		bool hasN = reads_pack[i].first;
		
		bool acceptRefRead = !hasN;
		if (referenceReadsMode == ReferenceReadsMode::Sparse)
			acceptRefRead &= ref_reads_accepter.ShouldAddToReference(current_read_id);

		if (acceptRefRead)
		{
			reference_reads.Add(*read);
			++id_in_reference;
		}

//		std::unordered_map<uint32_t, uint32_t> m;

/*		for (auto& kmer : read_kmers)
		{
			auto [localit, loc_end] = kmers.find(kmer);
			uint32_t kmer_card = 0;

			for (; localit != loc_end; ++localit)
			{
				++m[localit->second];
				++kmer_card;
			}

			if (acceptRefRead && kmer_card < maxKmerCount)
			{
				kmers.insert(kmer, id_in_reference - 1);
			}
		}*/

		m.clear();

		const uint32_t pf_prefix_offset = 2;
		const uint32_t pf_suffix_offset = 4;

		uint32_t read_kmers_size = static_cast<uint32_t>(read_kmers.size());
		for (uint32_t i = 0; i < read_kmers_size; ++i)
		{
			if(i + pf_prefix_offset < read_kmers_size)
				kmers.prefetch_prefix(read_kmers[i + pf_prefix_offset]);
			if(i + pf_suffix_offset < read_kmers_size)
				kmers.prefetch_suffix(read_kmers[i + pf_suffix_offset]);

			auto [localit, loc_end] = kmers.find(read_kmers[i]);
			uint32_t kmer_card = 0;

			for (; localit != loc_end; ++localit)
			{
				++m[localit->second];
				++kmer_card;
			}

			if (acceptRefRead && kmer_card < maxKmerCount)
			{
				kmers.insert(read_kmers[i], id_in_reference - 1);
			}

		}

		current_out_queue_elem.data.emplace_back();
		current_out_queue_elem.data.back().read_id = current_read_id;
		current_out_queue_elem.data.back().hasN = hasN;
		current_out_queue_elem.data.back().read = move(*read);

		neighbours.clear();
		neighbours.reserve(m.size());

		for (const auto& e : m)
			neighbours.emplace_back(e);

		auto min_iter = std::begin(neighbours);
		min_iter += std::min(max_candidates, (uint32_t) neighbours.size());

		std::partial_sort(std::begin(neighbours), min_iter, std::end(neighbours), [](const auto& e1, const auto& e2) {
			if (e1.second != e2.second)
				return e1.second > e2.second;
			return e1.first < e2.first;
		});

		if (neighbours.size() > max_candidates)
			neighbours.resize(max_candidates);

		for (const auto&  elem : neighbours)
			current_out_queue_elem.data.back().ref_reads.push_back(elem.first);

		++current_read_id;
	}

	//kmers.PrintMemoryUsage();
}

void CReadsSimilarityGraph::processReadsPackHiFi(const read_pack_t& reads_pack)
{
	auto accepted_kmers = getAcceptedKmers(reads_pack);

	std::vector<std::pair<uint32_t, std::pair<uint32_t, kmer_type>>> neighbours;

	std::vector<std::vector<kmer_type>> common_kmers;
	std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> m; //maps reads id to (number of commont k-mers with current read, id in commont k-mers container)

	m.max_load_factor(1.0);

	assert(reads_pack.size() == accepted_kmers.size());
	for(size_t i = 0 ; i < reads_pack.size() ; ++i)
	{
		auto& [read, read_kmers] = accepted_kmers[i];
		bool hasN = reads_pack[i].first;
		
		bool acceptRefRead = !hasN;
		if (referenceReadsMode == ReferenceReadsMode::Sparse)
			acceptRefRead &= ref_reads_accepter.ShouldAddToReference(current_read_id);

		if (acceptRefRead)
		{
			reference_reads.Add(*read);
			++id_in_reference;
		}
		
		m.clear();
		common_kmers.clear();

		const uint32_t pf_prefix_offset = 2;
		const uint32_t pf_suffix_offset = 4;

		uint32_t read_kmers_size = static_cast<uint32_t>(read_kmers.size());
		for (uint32_t i = 0; i < read_kmers_size; ++i)
		{
			if(i + pf_prefix_offset < read_kmers_size)
				kmers.prefetch_prefix(read_kmers[i + pf_prefix_offset]);
			if(i + pf_suffix_offset < read_kmers_size)
				kmers.prefetch_suffix(read_kmers[i + pf_suffix_offset]);

			auto [localit, loc_end] = kmers.find(read_kmers[i]);
			uint32_t kmer_card = 0;

			for (; localit != loc_end; ++localit)
			{
				auto& info = m[localit->second];
				++info.first;
				if (info.first == 1)
				{
					common_kmers.emplace_back();
					info.second = static_cast<uint32_t>(common_kmers.size()) - 1;
				}
				common_kmers[info.second].push_back(read_kmers[i]);

				++m[localit->second].first;
				++kmer_card;
			}

			if (acceptRefRead && kmer_card < maxKmerCount)
			{
				kmers.insert(read_kmers[i], id_in_reference - 1);
			}

		}

		current_out_queue_elem.data.emplace_back();
		current_out_queue_elem.data.back().read_id = current_read_id;
		current_out_queue_elem.data.back().hasN = hasN;
		current_out_queue_elem.data.back().read = move(*read);

		neighbours.clear();
		neighbours.reserve(m.size());

		for (const auto& e : m)
			neighbours.emplace_back(e);

		auto min_iter = std::begin(neighbours);
		min_iter += std::min(max_candidates, (uint32_t) neighbours.size());

		std::partial_sort(std::begin(neighbours), min_iter, std::end(neighbours), [](const auto& e1, const auto& e2) {
			if (e1.second.first != e2.second.first)
				return e1.second.first > e2.second.first;
			return e1.first < e2.first;
		});

		if (neighbours.size() > max_candidates)
			neighbours.resize(max_candidates);

		for (const auto& elem : neighbours)
		{
			current_out_queue_elem.data.back().ref_reads.push_back(elem.first);
			current_out_queue_elem.data.back().common_kmers.emplace_back(std::move(common_kmers[elem.second.second]));
		}

		++current_read_id;
	}

	//kmers.PrintMemoryUsage();
}

CReadsSimilarityGraph::CReadsSimilarityGraph(CParallelQueue<read_pack_t>& reads_queue,
	CParallelQueuePopWaiting<CCompressPack>& compress_queue,
	CReferenceReads& reference_reads,
	CReferenceGenome* reference_genome,
	const CKmerFilter& filteredKmers, uint32_t kmer_len, uint32_t max_candidates, uint32_t maxKmerCount, ReferenceReadsMode referenceReadsMode, CRefReadsAccepter& ref_reads_accepter,
	double refReadsFraction,
	int n_compression_threads,
	DataSource dataSource,
	double fill_factor_kmers_to_reads,
	bool verbose):
	compress_queue(compress_queue),
	reference_reads(reference_reads),
	filteredKmers(filteredKmers),
	kmer_len(kmer_len),
	max_candidates(max_candidates),
	maxKmerCount(maxKmerCount),
	kmers(kmer_len, static_cast<uint32_t>(filteredKmers.GetTotalKmers() * refReadsFraction), fill_factor_kmers_to_reads),
	referenceReadsMode(referenceReadsMode),
	ref_reads_accepter(ref_reads_accepter),
	n_compression_threads(n_compression_threads),
	dataSource(dataSource)
#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
	,
	internalThreads(std::make_unique<CReadsSimilarityGraphInternalThreads>(n_compression_threads + 1, kmer_len, filteredKmers, compress_queue))
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH
{
	processReferenceGenome(reference_genome);
	read_pack_t reads_pack;
	while (reads_queue.Pop(reads_pack))
	{
		if(dataSource == DataSource::PBHiFi)
			processReadsPackHiFi(reads_pack);
		else
			processReadsPack(reads_pack);

		current_out_queue_elem.id = current_out_elem_id++;
		compress_queue.Push(std::move(current_out_queue_elem));
	}

	if(verbose)
		kmers.PrintMemoryUsage();
	compress_queue.MarkCompleted();
#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
	internalThreads->NotifyFinish();
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH
}


#ifdef MEASURE_THREADS_TIMES
std::vector<double>& CReadsSimilarityGraph::GetTwTimes()
{
#ifdef USE_BETTER_PARALLELIZATION_IN_GRAPH
	return internalThreads->GetTwTimes();
#else
	return twTimes;
#endif // USE_BETTER_PARALLELIZATION_IN_GRAPH

}
#endif

CReadsSimilarityGraph::~CReadsSimilarityGraph() = default;
