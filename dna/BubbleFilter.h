#ifndef BUBBLEFilter_H_
#define BUBBLEFilter_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h" //e.g. for DefaultHash, and other serialization functions
#include "GlobalDna.h"

class BubbleSet
{
public:
	kmer_pair group;
	vector<Contig *> members;

	void free_members()
	{
		for(int i=0; i<members.size(); i++) delete members[i];
	}

	inline void add(Contig * v) //take
	{
		members.push_back(v);
	}

	inline void merge(BubbleSet * v) //copy and release
	{
		vector<Contig *> & mem = v->members;
		members.insert(members.end(), mem.begin(), mem.end());
		delete v;
	}

	friend ibinstream& operator<<(ibinstream& m, const BubbleSet& v)
	{
		m << v.group;
		m << v.members;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, BubbleSet& v)
	{
		m >> v.group;
		m >> v.members;
		return m;
	}
};

struct BubbleSet_hash
{
	size_t operator()(BubbleSet* x) const
	{
		kmer_pair & group = x->group;
		size_t seed = 0;
		hash_combine(seed, group.v1);
		hash_combine(seed, (group.v1 >> 32));
		hash_combine(seed, group.v2);
		hash_combine(seed, (group.v2 >> 32));
		return seed;
	}
};

struct BubbleSet_equal
{
public:
	bool operator()(BubbleSet* a,  BubbleSet* b) const
	{
		return a->group == b->group;
	}
};

class BubbleWorker
{
public:
	typedef hash_set<BubbleSet*, BubbleSet_hash, BubbleSet_equal> Groups;
	typedef Groups::iterator GroupIter;
	typedef vector<BubbleSet*> GroupVector;
	typedef vector<Contig *> ContigVector;

	KmerPairHash hash;
	Groups groups;
	GroupVector groupVec;
	ContigVector contigVec; //two parts: (1) dangling or long contigs (before shuffling), (2) non-filtered contigs (after shuffling)
	int edist_threshold;
	int len_threshold;

	BubbleWorker(int k, int edit_dist = 2, int length = 10000)
	{
		set_mer_length(k);
		edist_threshold = edit_dist;
		len_threshold = length;
	}

	~BubbleWorker()
	{
		for(int i=0; i < contigVec.size(); i++) delete contigVec[i];
	}

	void add_group(BubbleSet* group)
	{
		GroupIter itr = groups.find(group);
		if(itr != groups.end())
		{
			BubbleSet* cur = *itr;
			cur->merge(group);
		}
		else groups.insert(group);
	}

	void add_vertex(char* line)
	{
		Contig * v = new Contig;
		v->parse(line);
		if(v->in_neighbor != NULL_MER && v->out_neighbor != NULL_MER)
		{
			BubbleSet* grp = new BubbleSet;
			v->mark = grp->group.set(v->in_neighbor, v->out_neighbor);
			grp->add(v);
			add_group(grp);
		}
		else contigVec.push_back(v);
	}

	void sync_groups()
	{
		//distribute groups to sending buffers
		vector<GroupVector> _loaded_parts(_num_workers);
		for(GroupIter p = groups.begin(); p != groups.end(); p++)
		{
			BubbleSet* cur = *p;
			_loaded_parts[hash(cur->group)].push_back(cur);
		}
		groups.clear();
		//shuffle
		delete_after_all_to_all(_loaded_parts);
		//reset "groups"
		for(int i = 0; i < _num_workers; i++)
		{
			GroupVector & vec = _loaded_parts[i];
			for(int j = 0; j < vec.size(); j++) add_group(vec[j]);
		}
	};

	void groups_to_groupVec()
	{
		for(GroupIter p = groups.begin(); p != groups.end(); p++)  groupVec.push_back(*p);
		groups.clear();
	};

	void bubble_filter(BubbleSet* grp)
	{
		int size = grp->members.size();
		vector<bool> prune(size, false);
		for(int i=0; i<size; i++)
		{
			if(prune[i]) continue;
			Contig * a = grp->members[i];
			for(int j=i+1; j<size; j++)
			{
				Contig * b = grp->members[j];
				int editDis;
				if( a->mark != b->mark)
					editDis = (a->seq.length < b->seq.length) ? edist(a->seq.get_reverse(), b->seq) : edist(a->seq, b->seq.get_reverse());
				else
					editDis = edist(a->seq, b->seq);
				if( editDis <= edist_threshold)
				{
					if(a->freq > b->freq) prune[j] = true;
					else
					{
						prune[i] = true;
						break;
					}
				}
			}
		}
		//-----------------------
		for(int i=0; i<size; i++)
		{
			Contig * cur = grp->members[i];
			if(cur->seq.length <= len_threshold && prune[i]) delete cur;
			else contigVec.push_back(cur);
		}
	}

	//==============================

	void load_vertices(const char* inpath)
	{
		hdfsFS fs = getHdfsFS();
		hdfsFile in = getRHandle(inpath, fs);
		LineReader reader(fs, in);
		while (true)
		{
			reader.readLine();
			if (!reader.eof())
				add_vertex(reader.getLine());
			else
				break;
		}
		hdfsCloseFile(fs, in);
		hdfsDisconnect(fs);
	}

	//for machine read
	void dump_partition(const char* outpath)
	{
		hdfsFS fs = getHdfsFS();
		BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
		for(int i = 0 ; i < contigVec.size(); i++) contigVec[i]->dumpTo(writer);
		delete writer;
		hdfsDisconnect(fs);
	}

//        //for human read
//    	void dump_partition(const char* outpath)
//    	{
//        	hdfsFS fs = getHdfsFS();
//        	BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
//    		char buf[100];
//        	for(int i = 0 ; i < contigVec.size(); i++) {
//    			writer->check();
//    			writer->write(contigVec[i]->seq.toString().c_str());
//    			sprintf(buf, " %d\n", contigVec[i]->seq.length);
//    			writer->write(buf);
//        	}
//        	delete writer;
//        	hdfsDisconnect(fs);
//    	}

	//=======================================================
	// run the worker
	void run(const WorkerParams& params)
	{
		//check path + init
		if (_my_rank == MASTER_RANK)
		{
			if (dirCheck(params.input_path.c_str(), params.output_path.c_str(), _my_rank == MASTER_RANK, params.force_write) == -1)
				exit(-1);
		}
		init_timers();

		ResetTimer(WORKER_TIMER);
		vector<vector<string> >* arrangement;
		if (_my_rank == MASTER_RANK)
		{
			arrangement = params.native_dispatcher ? dispatchLocality(params.input_path.c_str()) : dispatchRan(params.input_path.c_str());

			masterScatter(*arrangement);
			vector<string>& assignedSplits = (*arrangement)[0];

			for (vector<string>::iterator it = assignedSplits.begin();
			        it != assignedSplits.end(); it++)
				load_vertices(it->c_str());
			delete arrangement;
		}
		else
		{
			vector<string> assignedSplits;
			slaveScatter(assignedSplits);

			for (vector<string>::iterator it = assignedSplits.begin();
			        it != assignedSplits.end(); it++)
				load_vertices(it->c_str());
		}
		StopTimer(WORKER_TIMER);
		PrintTimer("Load Time", WORKER_TIMER);

		ResetTimer(WORKER_TIMER);
		sync_groups();
		groups_to_groupVec();
		StopTimer(WORKER_TIMER);
		PrintTimer("Sync Time", WORKER_TIMER);

		ResetTimer(WORKER_TIMER);
		for(int i = 0; i < groupVec.size(); i++)
		{
			BubbleSet* cur = groupVec[i];
			bubble_filter(cur);
			delete cur;
		}
		groupVec.clear();
		StopTimer(WORKER_TIMER);
		PrintTimer("Bubble-Filter Time", WORKER_TIMER);

		ResetTimer(WORKER_TIMER);
		dump_partition(params.output_path.c_str());
		StopTimer(WORKER_TIMER);
		PrintTimer("Dump Time", WORKER_TIMER);
	}
};


void Bubble_Filter(string in_path, string out_path, int kmer, int editDis)
{
	WorkerParams bubble_p;
	bubble_p.input_path= in_path;
	bubble_p.output_path = out_path;
	bubble_p.force_write=true;
	bubble_p.native_dispatcher=false;
	//	init_workers();
	BubbleWorker bubble(kmer, editDis);  //second param means editDis
	bubble.run(bubble_p);
	//	worker_finalize();
}
#endif
