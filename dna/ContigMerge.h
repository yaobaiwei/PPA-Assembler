#ifndef CONTIGMERGE_H
#define CONTIGMERGE_H

#include "utils/type.h"
#include "basic/DNAPregel-dev.h" //e.g. for DefaultHash, and other serialization functions
#include "GlobalDna.h"

class ContigVertex
{
public:
	k_mer id;
	u8 type; //v1 = 1, v1-1 = 2
	neighbor_info neighbor1;
	neighbor_info neighbor2;
	u32 count1;
	u32 count2;

	k_mer parse(char* line) //line format: vid \t smaller_pred type(i.e. num_nbs) nb1 freq1 (nb2 freq2) ...
	{
		//return group_id
		char * pch;
		pch = strtok(line, "\t");
		id = strtoull(pch, NULL, 10);
		pch = strtok(NULL, " ");
		k_mer group = strtoull(pch, NULL, 10);
		pch = strtok(NULL, " ");
		type = atoi(pch);
		//------
		pch=strtok(NULL, " ");
		neighbor1.bitmap = atoi(pch); //highest 3 bits are 0, cannot be negative
		pch=strtok(NULL, " ");
		count1 = atoi(pch);
		//------
		if(type == 2)
		{
			pch=strtok(NULL, " ");
			neighbor2.bitmap = atoi(pch); //highest 3 bits are 0, cannot be negative
			pch=strtok(NULL, " ");
			count2 = atoi(pch);
		}
		return group;
	}

	void order_edges()
	{
		if(!neighbor1.dead_end()) if(!neighbor1.is_in()) neighbor1.reverse();
		if(!neighbor2.dead_end()) if(neighbor2.is_in()) neighbor2.reverse();
	}

	friend ibinstream& operator<<(ibinstream& m, const ContigVertex& v)
	{
		m << v.id;
		m << v.type;
		m << v.neighbor1;
		m << v.count1;
		if(v.type == 2)
		{
			m << v.neighbor2;
			m << v.count2;
		}
		return m;
	}

	friend obinstream& operator>>(obinstream& m, ContigVertex& v)
	{
		m >> v.id;
		m >> v.type;
		m >> v.neighbor1;
		m >> v.count1;
		if(v.type == 2)
		{
			m >> v.neighbor2;
			m >> v.count2;
		}
		return m;
	}
};

//=======================================

class ContigVSet
{
public:
	k_mer group;
	vector<ContigVertex *> members;

	void free_members()
	{
		for(int i=0; i<members.size(); i++) delete members[i];
	}

	inline void add(ContigVertex * v) //take
	{
		members.push_back(v);
	}

	inline void merge(ContigVSet * v) //copy and release
	{
		vector<ContigVertex *> & mem = v->members;
		members.insert(members.end(), mem.begin(), mem.end());
		delete v;
	}

	friend ibinstream& operator<<(ibinstream& m, const ContigVSet& v)
	{
		m << v.group;
		m << v.members;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, ContigVSet& v)
	{
		m >> v.group;
		m >> v.members;
		return m;
	}

	k_mer find_head_contig()
	{
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		ContigVertex * cur = members[member_pos[group]];
		if(cur->type==1)
			return cur->id;
		k_mer nb = get_neighbor(cur->id, cur->neighbor1);
		if(member_pos.find(nb) == member_pos.end())
			return cur->id;
		ContigVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			cur = members[member_pos[nb]];
			if(cur->type==1)
			{
				return cur->id;
			}
			nb = get_neighbor(cur->id, cur->neighbor1);
			if(nb == pre->id)
			{
				nb = get_neighbor(cur->id, cur->neighbor2);
			}
			if(member_pos.find(nb) == member_pos.end())
			{
				return cur->id;
			}
			pre = cur;
		}
		return NULL_MER;
	}

	k_mer find_tail_contig()
	{
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		ContigVertex * cur = members[member_pos[group]];
		k_mer nb;
		if(cur->type==1)
			nb =  get_neighbor(cur->id, cur->neighbor1);
		else
			nb = get_neighbor(cur->id, cur->neighbor2);
		if(member_pos.find(nb) == member_pos.end())
			return  cur->id;
		ContigVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			cur = members[member_pos[nb]];
			if(cur->type==1)
			{
				return  cur->id;
			}
			nb = get_neighbor(cur->id, cur->neighbor1);
			if(nb == pre->id)
			{
				nb = get_neighbor(cur->id, cur->neighbor2);
			}
			if(member_pos.find(nb) == member_pos.end())
			{
				return  cur->id;
			}
			pre = cur;
		}
		return NULL_MER;
	}

	Contig* get_contig()
	{
		//"id" of returned contig is meaningless, will be set later
		//remember to call free_members() after getting the contig (and then delete the current object itself)
		Contig * contig = new Contig;
		//================== A. create {member_ID -> member_pos} mapping ==================
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		//================== B. reorder members to "reordered", and adjust edge directions ==================
		vector<ContigVertex *> reordered(count);
		//--- 1. process the first vertex
		ContigVertex * cur = members[member_pos[group]]; //take the vertex whose ID equals group_ID
		reordered[0] = cur;
		if(cur->type == 1)
		{
			//v1 type, in_edge is the dead end
			//move old nb1 to new nb2 as out-edge
			cur->neighbor2 = cur->neighbor1;
			cur->count2 = cur->count1;
			//nb1 is the in-edge, i.e., dead end
			cur->neighbor1.bitmap = DEAD_END;
			cur->count1 = 0;
		}
		else
		{
			//pick a v{n-m}/dead_end neighbor as in-neighbor
			k_mer nb_in = get_neighbor(cur->id, cur->neighbor1);
			if(member_pos.find(nb_in) != member_pos.end())
			{
				//neighbor2 should be in-neighbor, swap two neighbors
				neighbor_info tmp = cur->neighbor1;
				cur->neighbor1 = cur->neighbor2;
				cur->neighbor2 = tmp;
				//---
				u32 tmp1 = cur->count1;
				cur->count1 = cur->count2;
				cur->count2 = tmp1;
			}
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//If the contig is a loop, the neighbor1 of the first node also can be one of the members
			nb_in = get_neighbor(cur->id, cur->neighbor1);
			if(member_pos.find(nb_in) != member_pos.end())
			{
				cur->neighbor1.bitmap =  DEAD_END;
				cur->count1 = 0;
			}
		}
		cur->order_edges();
		//--- 2. process other vertices
		ContigVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			k_mer curID = get_neighbor(pre->id, pre->neighbor2);
			cur = members[member_pos[curID]];
			reordered[i] = cur;
			//--- to make "pre" as "neighbor1"
			k_mer nb1_ID = get_neighbor(curID, cur->neighbor1);
			if(nb1_ID != pre->id)
			{
				neighbor_info tmp = cur->neighbor1;
				cur->neighbor1 = cur->neighbor2;
				cur->neighbor2 = tmp;
				//---
				u32 tmp1 = cur->count1;
				cur->count1 = cur->count2;
				cur->count2 = tmp1;
			}
			//------
			cur->order_edges();
			pre = cur;
		}
		//--- 3. handle the case when last vertex is a tip
		if(count > 1)
		{
			if(reordered.back()->type == 1)
			{
				cur->neighbor2.bitmap = DEAD_END;
				cur->count2 = 0;
			}
			else
			{
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//Also should consider the loop-contig case for the last node's neighbor
				k_mer nb_out = get_neighbor(cur->id, cur->neighbor2);
				if(member_pos.find(nb_out) != member_pos.end())
				{
					cur->neighbor2.bitmap =  DEAD_END;
					cur->count2 = 0;
				}
			}
		}

		//================== C. construct contig ==================
		//--- 1. set in_neighbor, in_pol, in_count
		neighbor_info & left_end = reordered.front()->neighbor1;
		if(!left_end.dead_end())
		{
			contig->in_neighbor = get_neighbor(reordered.front()->id, left_end);
			contig->in_pol = left_end.left_isH();
			contig-> in_count = reordered.front()->count1;
		}
		else
		{
			contig->in_neighbor = NULL_MER;
			contig->in_pol = false;
			contig-> in_count = 0;
		}
		//--- 2. set out_neighbor, out_pol, out_count
		neighbor_info & right_end = reordered.back()->neighbor2;
		if(!right_end.dead_end())
		{
			contig->out_neighbor = get_neighbor(reordered.back()->id, right_end);
			contig->out_pol = right_end.right_isH();
			contig->out_count = reordered.back()->count2;
		}
		else
		{
			contig->out_neighbor = NULL_MER;
			contig->out_pol = false;
			contig->out_count = 0;
		}
		//--- 3. set seq
		ATGC_bitmap & seq = contig->seq;
		bool pol; //whether first vertex is H?
		if(!left_end.dead_end()) pol = left_end.right_isH(); //look at in-edge
		else pol = reordered.front()->neighbor2.left_isH(); //look at out-edge
		if(pol) seq.init(getRC(reordered.front()->id)); //first vertex is H
		else seq.init(reordered.front()->id); //first vertex is L
		u32 min_freq = UINT_MAX;
		for(int i=1; i<count; i++)
		{
			ContigVertex * v = reordered[i];
			if(v->count1 < min_freq) min_freq = v->count1;
			pol = v->neighbor1.right_isH(); //look at in-edge
			if(pol)
			{
				k_mer to_append = v->id >> (mer_length*2 -2);
				seq.append(to_append ^ 3ull); //complement of first tag
			}
			else seq.append(v->id & 3ull); //last tag
		}
		contig->freq = min_freq; //may be UINT_MAX if it's a singleton contig
		return contig;
	}
};

//ContigVSet Key Order  =====================
//for defining hash_set<ContigVSet>
struct ContigVSet_hash
{
	size_t operator()(ContigVSet* x) const
	{
		return x->group;
	}
};

struct ContigVSet_equal
{
public:
	bool operator()(ContigVSet* a,  ContigVSet* b) const
	{
		return a->group == b->group;
	}
};

bool ContigVSort( ContigVSet* m1,  ContigVSet* m2)
{
	return m1->group  < m2->group;
}

//=======================================
static int tip_threshold;

class ContigWorker
{
public:
	typedef hash_set<ContigVSet*, ContigVSet_hash, ContigVSet_equal> Groups;
	typedef Groups::iterator GroupIter;
	typedef vector<ContigVSet*> GroupVector;
	typedef vector<Contig*> ContigVector;

	DefaultHash<k_mer> hash;
	Groups groups;
	GroupVector groupVec;
	ContigVector contigs;

	ContigWorker(int k, int tip_length = 0)
	{
		set_mer_length(k);
		tip_threshold = tip_length;
	}

	~ContigWorker()
	{
		for(int i = 0; i < contigs.size(); i++) delete contigs[i];
	}

	void add_group(ContigVSet* group)
	{
		GroupIter itr = groups.find(group);
		if(itr != groups.end())
		{
			ContigVSet* cur = *itr;
			cur->merge(group);
		}
		else groups.insert(group);
	}

	void add_vertex(char* line)
	{
		ContigVertex * v = new ContigVertex;
		ContigVSet* grp = new ContigVSet;
		grp->group = v->parse(line);
		grp->add(v);
		add_group(grp);
	}

	void sync_groups()
	{
		//distribute groups to sending buffers
		vector<GroupVector> _loaded_parts(_num_workers);
		for(GroupIter p = groups.begin(); p != groups.end(); p++)
		{
			ContigVSet* cur = *p;
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
		//cout<<"Worker "<<_my_rank<<": \""<<inpath<<"\" loaded"<<endl;//DEBUG !!!!!!!!!!
	}

	//* //for machine read
	void dump_partition(const char* outpath)
	{
		//id \t in_nb out_nb freq seq
		//in_nb/out_nb = (nb_id, pol, count) or (NULL_MER, 0, 0)
		//seq = (length num_bytes byte1 byte2 ...)
		hdfsFS fs = getHdfsFS();
		BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
		for(int i = 0; i < contigs.size(); i++)
		{
			Contig * cur = contigs[i];
			cur->dumpTo(writer);
		}
		delete writer;
		hdfsDisconnect(fs);
	}
	//*/


//    //for human read
//	void dump_partition(const char* outpath)
//	{//one contig per line: contig length
//		hdfsFS fs = getHdfsFS();
//		BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
//		char buf[100];
//		for(int i = 0; i < contigs.size(); i++) {
//			writer->check();
//			writer->write(contigs[i]->seq.toString().c_str());
//			sprintf(buf, " %d\n", contigs[i]->seq.length);
//			writer->write(buf);
//		}
//		delete writer;
//		hdfsDisconnect(fs);
//	}


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
				load_vertices(it->c_str());
			delete arrangement;
		}
		else
		{
			vector<string> assignedSplits;
			slaveScatter(assignedSplits);
			//reading assigned splits (map)
			for (vector<string>::iterator it = assignedSplits.begin();
			        it != assignedSplits.end(); it++)
				load_vertices(it->c_str());
		}
		StopTimer(WORKER_TIMER);
		PrintTimer("Load Time", WORKER_TIMER);

		//send vertex groups according to hash_id (reduce)
		ResetTimer(WORKER_TIMER);
		sync_groups();
		groups_to_groupVec();
		StopTimer(WORKER_TIMER);
		PrintTimer("Sync Time", WORKER_TIMER);

		//merge each group into a contig
		ResetTimer(WORKER_TIMER);

#ifdef SV_USED
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//add the following function in the case that we used SV-alg in previous step,
		//in that way, the group is not the head of this contig, we should find the head
		for(int i = 0; i < groupVec.size(); i++)
		{
			ContigVSet* cur = groupVec[i];
			k_mer head = cur->find_head_contig();
			k_mer tail = cur->find_tail_contig();
			cur->group = (head > tail) ? tail : head;
		}
#endif
		sort(groupVec.begin(), groupVec.end(), ContigVSort);

		for(int i = 0; i < groupVec.size(); i++)
		{
			ContigVSet* cur = groupVec[i];
			Contig * tmp = cur->get_contig();
			cur->free_members();
			delete cur;
			if( (tmp->seq.length <= tip_threshold) && ((tmp->in_neighbor == NULL_MER) || (tmp->out_neighbor == NULL_MER)) )
			{
				delete tmp;
				continue;
			}
			//set id
			tmp->id = _my_rank;
			tmp->id <<= 32;
			tmp->id |= (NULL_MER | i);
			contigs.push_back(tmp);
		}
		StopTimer(WORKER_TIMER);
		PrintTimer("Merge Time", WORKER_TIMER);

		//dump De Bruijn graph
		ResetTimer(WORKER_TIMER);
		dump_partition(params.output_path.c_str());
		StopTimer(WORKER_TIMER);
		PrintTimer("Dump Time", WORKER_TIMER);
	}
};


void Contig_Merge(string in_path, string out_path, int kmer)
{
	WorkerParams contig_p;
	contig_p.input_path= in_path;
	contig_p.output_path = out_path;
	contig_p.force_write=true;
	contig_p.native_dispatcher=false;
	//	init_workers();
	ContigWorker contig(kmer);
	contig.run(contig_p);
	//	worker_finalize();
}

#endif
