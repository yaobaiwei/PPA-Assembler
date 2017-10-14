#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <vector>
#include <string>

#include "utils/time.h"
#include "utils/global.h"
#include "utils/communication.h"
#include "utils/serialization.h"
#include "utils/hdfs_core.h"
#include "utils/combiner.h"
#include "utils/aggregator.h"
#include "utils/type.h"
#include "basic/Vertex.h"
#include "GlobalDna.h"
using namespace std;

k_mer ALL_A, ALL_C, ALL_G, ALL_T;

//====================================
class KPlus_mer
{
public:
	k_mer id;
	u8 count;

	KPlus_mer()
	{
		id = 0;
		count = 0;
	}

	//operations on vid
	void append2id(char c)
	{
		id <<= 2;
		if (c == 'T')
			id |= 3ull;
		else if (c == 'G')
			id |= 2ull;
		else if (c == 'C')
			id |= 1ull;
	}

	void id2canonical(char * buf)
	{
		char rev[mer_length + 1];
		for (int i = 0; i < mer_length + 1; i++)
		{
			char tmp = buf[mer_length - i];
			if (tmp == 'A')
				rev[i] = 'T';
			else if (tmp == 'T')
				rev[i] = 'A';
			else if (tmp == 'G')
				rev[i] = 'C';
			else if (tmp == 'C')
				rev[i] = 'G';
		}
		if (strncmp(buf, rev, mer_length + 1) < 0)
		{
			for (int i = 0; i < mer_length + 1; i++)
				append2id(buf[i]);
		}
		else
		{
			for (int i = 0; i < mer_length + 1; i++)
				append2id(rev[i]);
		}
	}

	inline k_mer get_left_kmer()
	{
		return (id >> 2);
	}

	inline k_mer get_right_kmer()
	{
		return (id & kick);//kick out highest two bits
	}

	inline u32 get_leftmost()
	{
		return (u32)(id >> 2 * mer_length);
	}

	inline u32 get_rightmost()
	{
		return  (u32)(id & 3ull);
	}

	friend ibinstream& operator<<(ibinstream& m, const KPlus_mer& v)
	{
		m << v.id;
		m << v.count;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, KPlus_mer& v)
	{
		m >> v.id;
		m >> v.count;
		return m;
	}
};

//k+1 mer Key Order  =====================
//for defining hash_set<k+1 mer>
struct vid_hash_kplus
{
	size_t operator()(KPlus_mer* x) const
	{
		return x->id;
	}
};

struct vid_equal_kplus
{
public:
	bool operator()(KPlus_mer* a, KPlus_mer* b) const
	{
		if (a->id == b->id)
			return true;
		else
			return false;
	}
};
//====================================

class DNAVertex
{
public:
	k_mer id;
	u32 bitmap; //four bytes, each for LL, LH, HL, and HH
	vector<u8> freqs;

	DNAVertex()
	{
		id = 0;
		bitmap = 0;
	}

	void set_edgeBit(u32 tag, bool is_in, int shift)
	{
		tag = bit_pos[tag];
		if(is_in)
			tag <<= 4;
		bitmap |= (tag << shift);
	}

	bool vid2canonical(k_mer vid)
	{
		k_mer rc = getRC(vid);

		if(rc > vid)
		{
			id = vid;
			return false;
		}
		else
		{
			id = rc;
			return true;
		}
	}

	friend ibinstream& operator<<(ibinstream& m, const DNAVertex& v)
	{
		m << v.id;
		m << v.bitmap;
		m << v.freqs;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, DNAVertex& v)
	{
		m >> v.id;
		m >> v.bitmap;
		m >> v.freqs;
		return m;
	}
};

//DNA Key Order  =====================
//for defining hash_set<DNAVertex>
struct vid_hash
{
	size_t operator()(DNAVertex* x) const
	{
		return x->id;
	}
};

struct vid_equal
{
public:
	bool operator()(DNAVertex* a,  DNAVertex* b) const
	{
		if(a->id == b->id)
			return true;
		else
			return false;
	}
};
//====================================

class DeBruijn
{
public:
	typedef hash_set<DNAVertex*, vid_hash, vid_equal> VertexContainer;
	typedef VertexContainer::iterator VertexIter;
	typedef vector<DNAVertex*> VertexVector;

	typedef hash_set<KPlus_mer*, vid_hash_kplus, vid_equal_kplus> KPlusContainer;
	typedef KPlusContainer::iterator KPlusIter;
	typedef vector<KPlus_mer*> KPlusVector;

	int freq_threshold;
	DefaultHash<k_mer> hash;
	VertexContainer vertexes;
	KPlusContainer kplus_mers;

	DeBruijn(int k, int freq)
	{
		set_mer_length(k);
		get_loop_kplus();
		freq_threshold = freq;
	}

	void freeVContainer()
	{
		//note: cannot iterate through the set "vertexes", as deletion will make v->id invalid
		VertexVector to_del;
		for (VertexIter p = vertexes.begin(); p != vertexes.end(); p++)
			to_del.push_back(*p);
		vertexes.clear();
		for (int i = 0; i<to_del.size(); i++)
			delete to_del[i];
	}

	~DeBruijn()
	{
		freeVContainer();
	}

	void get_loop_kplus()
	{
		ALL_A = 0;
		ALL_C = 0x5555555555555555 >> (62 - 2 * mer_length);
		ALL_G = 0xAAAAAAAAAAAAAAAA >> (62 - 2 * mer_length);
		ALL_T = 0xFFFFFFFFFFFFFFFF >> (62 - 2 * mer_length);
	}

	bool is_loop_kplus(KPlus_mer* kplus)
	{
		if(kplus->id == ALL_A || kplus->id == ALL_C || kplus->id == ALL_G || kplus->id == ALL_T )
			return true;
		return false;
	}

	void split(char * str, vector<char *> & out) //split by 'N'
	{
		char * cur = str;
		for (int i = 0; str[i] != '\0'; i++)
		{
			if (str[i] == 'N')
			{
				str[i] = '\0';
				out.push_back(cur);
				cur = str + i + 1;
			}
		}
		out.push_back(cur);
	}

	//==============================

	void add_kplus_mer(KPlus_mer* kplus)
	{
		if(is_loop_kplus(kplus))
		{
			delete kplus;
			return;
		}
		KPlusIter itr = kplus_mers.find(kplus);
		if (itr != kplus_mers.end())
		{
			KPlus_mer * v = *itr;
			v->count += kplus->count;
			delete kplus;
		}
		else
			kplus_mers.insert(kplus);
	}

	void add_kplus_mers(char* line)
	{
		vector<char *> reads;
		split(line, reads);
		for (int i = 0; i<reads.size(); i++)
		{
			char * read = reads[i];
			int length = strlen(read);
			if(length < mer_length + 1)
				continue;
			int count = length - mer_length;
			for(int j=0; j<count; j++)
			{
				char * kmer = read + j;
				KPlus_mer * kplus = new KPlus_mer;
				kplus->id2canonical(kmer);
				kplus->count = 1;
				add_kplus_mer(kplus);
			}
		}
	}

	//==============================

	void add_vertex(DNAVertex* vertex)
	{
		VertexIter itr = vertexes.find(vertex);
		if(itr != vertexes.end())
		{
			DNAVertex * v = *itr;
			//------
			vector<u32> my_count, new_count, merged;
			int my_pos = 0, new_pos = 0;
			parse_vints(my_count, v->freqs);
			parse_vints(new_count, vertex->freqs);
			//------
			bool v1_pol[4] = {false,false,true,true};
			bool v2_pol[4] = {false, true, false, true};

			for(int j = 0; j < 4; j++)
			{
				//LL LH HL HH
				int shift = getShift(v1_pol[j], v2_pol[j]);
				//------
				for(int i=0; i<8; i++)
				{
					u32 bit = ATGC_bits[i];
					bool me_exists = (v->bitmap >> shift) & bit;
					bool new_exists = (vertex->bitmap >> shift) & bit;

					if(me_exists && new_exists)
					{
						merged.push_back(my_count[my_pos++] + new_count[new_pos++]);
					}
					else if(me_exists)
					{
						merged.push_back(my_count[my_pos++]);
					}
					else if(new_exists)
					{
						v->bitmap |= (bit << shift);
						merged.push_back(new_count[new_pos++]);
					}
					else
					{}
				}
			}
			//------
			vector<u8>().swap(v->freqs);
			for(int i=0; i<merged.size(); i++)
			{
				to_vint(merged[i]);
				append_vint(v->freqs);
			}
			delete vertex;
		}
		else
			vertexes.insert(vertex);
	}

	void add_vertices(KPlusVector & vec)
	{
		for (int i = 0; i<vec.size(); i++)
		{
			KPlus_mer* kplus = vec[i];
			if (kplus->count < freq_threshold)
			{
				delete kplus;
				continue; //filter out low-freq k+1 mers
			}
			//------
			DNAVertex *v1 = new DNAVertex;
			bool v1_pol = v1->vid2canonical(kplus->get_left_kmer());
			DNAVertex *v2 = new DNAVertex;
			bool v2_pol = v2->vid2canonical(kplus->get_right_kmer());
			int shift = getShift(v1_pol, v2_pol);
			//------
			to_vint(kplus->count);
			v1->set_edgeBit(kplus->get_rightmost(), false, shift);
			append_vint(v1->freqs);
			v2->set_edgeBit(kplus->get_leftmost(), true, shift);
			append_vint(v2->freqs);
			//------
			add_vertex(v1);
			add_vertex(v2);
			delete kplus;
		}
	}

	void reduce_kplus_mers()
	{
		vector<KPlusVector> _loaded_parts(_num_workers);
		for(KPlusIter p = kplus_mers.begin(); p != kplus_mers.end(); p++)
		{
			KPlus_mer* kplus = *p;
			_loaded_parts[hash(kplus->id)].push_back(kplus);
		}
		//------
		kplus_mers.clear();
		KPlusContainer().swap(kplus_mers);
		//------
		delete_after_all_to_all(_loaded_parts);
		for(int i = 0; i < _num_workers; i++)
		{
			KPlusVector & vec = _loaded_parts[i];
			for (int j = 0; j < vec.size(); j++)
				add_kplus_mer(vec[j]);
		}
		vector<KPlusVector>().swap(_loaded_parts);
		//------
		KPlusVector kplus_vec;
		for (KPlusIter p = kplus_mers.begin(); p != kplus_mers.end(); p++)
			kplus_vec.push_back(*p);
		kplus_mers.clear();
		KPlusContainer().swap(kplus_mers);
		//------
		add_vertices(kplus_vec);
	};

	//==================================
	void sync_graph()
	{
		vector<VertexVector> _loaded_parts(_num_workers);
		for (VertexIter p = vertexes.begin(); p != vertexes.end(); p++)
		{
			DNAVertex* v = *p;
			_loaded_parts[hash(v->id)].push_back(v);
		}
		vertexes.clear();
		VertexContainer().swap(vertexes);
		delete_after_all_to_all(_loaded_parts);
		for (int i = 0; i < _num_workers; i++)
		{
			VertexVector & vec = _loaded_parts[i];
			for(int j = 0; j < vec.size(); j++)
				add_vertex(vec[j]);
		}
		_loaded_parts.clear();
	};

	void load_graph(const char* inpath)
	{
		hdfsFS fs = getHdfsFS();
		hdfsFile in = getRHandle(inpath, fs);
		LineReader reader(fs, in);
		while (true)
		{
			reader.readLine();
			if (!reader.eof())
				add_kplus_mers(reader.getLine());
			else
				break;
		}
		hdfsCloseFile(fs, in);
		hdfsDisconnect(fs);
		//cout<<"Worker "<<_my_rank<<": \""<<inpath<<"\" loaded"<<endl;//DEBUG !!!!!!!!!!
	}

	//==============================
	void dump_partition(const char* outpath)
	{
		hdfsFS fs = getHdfsFS();
		BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
		char buf[100];
		for (VertexIter it = vertexes.begin(); it != vertexes.end(); it++)
		{
			writer->check();
			//writer->write((*it)->toString().c_str());
			// Adding....
			DNAVertex * v = (*it);
			vector<u8> & vec = v->freqs;
			sprintf(buf, "%llu\t%u %d", v->id, v->bitmap, vec.size());
			writer->write(buf);
			for(int i=0; i<vec.size(); i++)
			{
				sprintf(buf, " %hhu", vec[i]);
				writer->write(buf);
			}
			//.................
			writer->write("\n");
		}
		delete writer;
		hdfsDisconnect(fs);
	}

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

		//dispatch splits
		ResetTimer(WORKER_TIMER);
		vector<vector<string> >* arrangement;
		if (_my_rank == MASTER_RANK)
		{
			arrangement = params.native_dispatcher ? dispatchLocality(params.input_path.c_str()) : dispatchRan(params.input_path.c_str());
			//reportAssignment(arrangement);//DEBUG !!!!!!!!!!
			masterScatter(*arrangement);
			vector<string>& assignedSplits = (*arrangement)[0];
			//reading assigned splits (map)
			for (vector<string>::iterator it = assignedSplits.begin();
			        it != assignedSplits.end(); it++)
				load_graph(it->c_str());
			delete arrangement;
		}
		else
		{
			vector<string> assignedSplits;
			slaveScatter(assignedSplits);
			//reading assigned splits (map)
			for (vector<string>::iterator it = assignedSplits.begin();
			        it != assignedSplits.end(); it++)
				load_graph(it->c_str());
		}
		StopTimer(WORKER_TIMER);
		PrintTimer("Load Time", WORKER_TIMER);

		//send vertices according to hash_id (reduce)
		ResetTimer(WORKER_TIMER);
		reduce_kplus_mers();
		sync_graph();

		//sync_graph();
		StopTimer(WORKER_TIMER);
		PrintTimer("Sync Time", WORKER_TIMER);

		//dump De Bruijn graph
		ResetTimer(WORKER_TIMER);
		dump_partition(params.output_path.c_str());
		StopTimer(WORKER_TIMER);
		PrintTimer("Dump Time", WORKER_TIMER);
	}
};


void DeBruijn_Build(string in_path, string outpath, int kmer, int freq_t)
{
	WorkerParams deBruijn_p;
	deBruijn_p.input_path= in_path;
	deBruijn_p.output_path= outpath;
	deBruijn_p.force_write=true;
	deBruijn_p.native_dispatcher=false;
	//	init_workers();
	DeBruijn deBruijn(kmer, freq_t);
	deBruijn.run(deBruijn_p);
	//	worker_finalize();
}

#endif
