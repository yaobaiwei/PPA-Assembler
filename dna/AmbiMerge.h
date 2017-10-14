#ifndef AMBIMERGE_H_
#define AMBIMERGE_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h" //e.g. for DefaultHash, and other serialization functions
#include "GlobalDna.h"
#include "DNAMessageBuffer.h"


class AmbiVertex
{
public:
	k_mer id;
	u8 type;
	vector<AmbiNB> ambi_nbs;
	vector<ContigNB> contig_nbs;

	k_mer parse(char* line)
	{
		//return group_id
		char * pch;
		pch = strtok(line, "\t");
		id = strtoull(pch, NULL, 10);
		pch = strtok(NULL, " ");
		k_mer group = strtoull(pch, NULL, 10);
		pch = strtok(NULL, " ");
		type = atoi(pch);
		pch=strtok(NULL, " ");
		int size = atoi(pch);
		for(int i = 0; i < size; i++)
		{
			pch=strtok(NULL, " ");
			AmbiNB ambiNB;
			ambiNB.nb_info.bitmap = (u8)atoi(pch);
			pch=strtok(NULL, " ");
			ambiNB.count = atoi(pch);
			ambi_nbs.push_back(ambiNB);
		}
		pch=strtok(NULL, " ");
		size = atoi(pch);
		for(int i = 0; i< size; i++)
		{
			pch=strtok(NULL, " ");
			ContigNB contigNB;
			contigNB.nid = strtoull(pch, NULL, 10);
			pch=strtok(NULL, " ");
			contigNB.ninfo.bitmap = (u8)atoi(pch);
			pch=strtok(NULL, " ");
			contigNB.contigID = strtoull(pch, NULL, 10);
			pch=strtok(NULL, " ");
			contigNB.length = atoi(pch);
			contig_nbs.push_back(contigNB);
		}
		return group;
	}

	friend ibinstream& operator<<(ibinstream& m, const AmbiVertex& v)
	{
		m << v.id;
		m << v.type;
		m << v.ambi_nbs;
		m << v.contig_nbs;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, AmbiVertex& v)
	{
		m >> v.id;
		m >> v.type;
		m >> v.ambi_nbs;
		m >> v.contig_nbs;
		return m;
	}
};

class AmbiContig
{
public:
	ATGC_bitmap seq;

	void dumpTo(BufferedWriter* writer, int line_num)
	{
		char buf[50000];
		writer->check();
		sprintf(buf, "%s-%d-%d\n", ">contig", _my_rank,line_num+1);
		writer->write(buf);
		sprintf(buf, "%s\n", seq.toString().c_str());
		writer->write(buf);
	}

	friend ibinstream& operator<<(ibinstream& m, const AmbiContig& v)
	{
		m << v.seq;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, AmbiContig& v)
	{
		m >> v.seq;
		return m;
	}
};

//=======================================
struct vtuple
{
	k_mer vid;
	bool isreverse;
	bool iscontig;

	vtuple() {}

	vtuple(k_mer id, bool is_reverse, bool is_contig)
	{
		vid = id;
		isreverse = is_reverse;
		iscontig = is_contig;
	}
};

ibinstream& operator<<(ibinstream& m, const vtuple& v)
{
	m << v.vid;
	m << v.isreverse;
	m << v.iscontig;
	return m;
}

obinstream& operator>>(obinstream& m, vtuple& v)
{
	m >> v.vid;
	m >> v.isreverse;
	m >> v.iscontig;
	return m;
}

struct MessageValue
{
	k_mer id;
	ATGC_bitmap sequence;
};

ibinstream & operator<<(ibinstream & m, const MessageValue & v)
{
	m << v.id;
	m << v.sequence;
	return m;
}

obinstream & operator>>(obinstream & m,  MessageValue & v)
{
	m >> v.id;
	m >> v.sequence;
	return m;
}

class AmbiVSet
{
public:
	k_mer id;
	vector<AmbiVertex *> members;
	vector<vtuple> order;

	void free_members()
	{
		for(int i=0; i<members.size(); i++) delete members[i];
	}

	inline void add(AmbiVertex * v) //take
	{
		members.push_back(v);
	}

	inline void merge(AmbiVSet * v) //copy and release
	{
		vector<AmbiVertex *> & mem = v->members;
		members.insert(members.end(), mem.begin(), mem.end());
		delete v;
	}

	friend ibinstream& operator<<(ibinstream& m, const AmbiVSet& v)
	{
		m << v.id;
		m << v.members;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, AmbiVSet& v)
	{
		m >> v.id;
		m >> v.members;
		return m;
	}

	k_mer find_head_ambicontig()
	{
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		AmbiVertex * cur = members[member_pos[id]];
		if(cur->type== V_1)
			return cur->id;
		k_mer next;
		if(cur->ambi_nbs.size() == 2)
		{
			next = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
		}
		else
		{
			next = cur->contig_nbs[0].nid;
		}
		if(member_pos.find(next) == member_pos.end())
			return cur->id;
		AmbiVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			cur = members[member_pos[next]];
			if(cur->type == V_1)
			{
				return cur->id;
			}
			if(cur->ambi_nbs.size() == 2)
			{
				next = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				if(next == pre->id)
				{
					next =  get_neighbor(cur->id, cur->ambi_nbs[1].nb_info);
				}
			}
			else if (cur->contig_nbs.size() == 2)
			{
				next = cur->contig_nbs[0].nid;
				if(next == pre->id)
				{
					next =  cur->contig_nbs[1].nid;
				}
			}
			else
			{
				next  = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				if(next == pre->id)
				{
					next  = cur->contig_nbs[0].nid;
				}
			}
			if(member_pos.find(next) == member_pos.end())
			{
				return cur->id;
			}
			pre = cur;
		}
		return NULL_MER;
	}

	k_mer find_tail_ambicontig()
	{
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		AmbiVertex * cur = members[member_pos[id]];
		k_mer next;
		if(cur->type== V_1)
		{
			if(cur->ambi_nbs.size())
			{
				neighbor_info & ninfo = cur->ambi_nbs[0].nb_info;
				next = get_neighbor(cur->id, ninfo);
			}
			else
			{
				next = cur->contig_nbs[0].nid;
			}
		}
		else
		{
			if(cur->ambi_nbs.size() == 2)
			{
				next = get_neighbor(cur->id, cur->ambi_nbs[1].nb_info);
			}
			else if (cur->contig_nbs.size() == 2)
			{
				next = cur->contig_nbs[1].nid;
			}
			else
			{
				next = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
			}
		}
		if(member_pos.find(next) == member_pos.end())
			return cur->id;
		AmbiVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			cur = members[member_pos[next]];
			if(cur->type == V_1)
			{
				return cur->id;
			}
			if(cur->ambi_nbs.size() == 2)
			{
				next = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				if(next == pre->id)
				{
					next =  get_neighbor(cur->id, cur->ambi_nbs[1].nb_info);
				}
			}
			else if (cur->contig_nbs.size() == 2)
			{
				next = cur->contig_nbs[0].nid;
				if(next == pre->id)
				{
					next =  cur->contig_nbs[1].nid;
				}
			}
			else
			{
				next  = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				if(next == pre->id)
				{
					next  = cur->contig_nbs[0].nid;
				}
			}
			if(member_pos.find(next) == member_pos.end())
			{
				return cur->id;
			}
			pre = cur;
		}
		return NULL_MER;
	}

	void  reorder_ambiVSet(k_mer head = NULL_MER)
	{
		hash_map<k_mer, int> member_pos;
		int count = members.size();
		for(int i = 0; i < count; i++) member_pos[members[i]->id] = i;
		//================== B. reorder members to "reordered" ==================
		vector<AmbiVertex *> reordered(count);
		//--- 1. process the first vertex
		head = (head == NULL_MER) ? id : head;
		AmbiVertex * cur = members[member_pos[head]]; //take the vertex whose ID equals group_ID
		reordered[0] = cur;
		k_mer next;
		if(cur->type == V_1)
		{
			// has only one ambi_nb
			if(cur->ambi_nbs.size())
			{
				neighbor_info & ninfo = cur->ambi_nbs[0].nb_info;
				k_mer nb = get_neighbor(cur->id, ninfo);
				if(ninfo.is_in())
					ninfo.reverse();
				order.push_back(vtuple(cur->id, ninfo.left_isH(), false));
				if(member_pos.find(nb) != member_pos.end())
				{
					next = nb;
				}
				else
				{
					next = NULL_MER;
				}
			}
			else
			{
				//has only one contig_nb
				k_mer nb = cur->contig_nbs[0].nid;
				if(member_pos.find(nb) != member_pos.end())
					next = nb;
				else next = NULL_MER;
				bool has_reverse = false;
				//we should reverse the direction if the contig is not in
				if(cur->contig_nbs[0].ninfo.is_in())
				{
					cur->contig_nbs[0].ninfo.reverse();
					has_reverse = true;
				}
				order.push_back(vtuple(cur->id, cur->contig_nbs[0].ninfo.left_isH(), false));
				order.push_back(vtuple(cur->contig_nbs[0].contigID, has_reverse, true));
			}
		}
		else
		{
			if(cur->ambi_nbs.size() == 2)
			{
				k_mer nb1 = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				k_mer nb2 = get_neighbor(cur->id, cur->ambi_nbs[1].nb_info);
				int index;
				if(member_pos.find(nb1) != member_pos.end())
				{
					next = nb1;
					index = 0;
				}
				else
				{
					next = nb2;
					index = 1;
				}
				if(cur->ambi_nbs[index].nb_info.is_in())
				{
					cur->ambi_nbs[index].nb_info.reverse();
				}
				order.push_back(vtuple(cur->id, cur->ambi_nbs[index].nb_info.left_isH(), false));
			}
			else if (cur->contig_nbs.size() == 2)
			{
				k_mer nb1 = cur->contig_nbs[0].nid;
				k_mer nb2 = cur->contig_nbs[1].nid;
				int left, right;
				if(member_pos.find(nb1) != member_pos.end())
				{
					next = nb1;
					left = 1;
					right = 0;
				}
				else if(nb1 == NULL_MER)
				{
					if(member_pos.find(nb2) != member_pos.end())
					{
						next = nb2;
						left = 0;
						right = 1;
					}
					else
					{
						next = NULL_MER;
						left = 1;
						right = 0;
					}
				}
				else
				{
					next = nb2;  //nb2 also can be NULL_MER in here
					left = 0;
					right = 1;
				}
				bool has_reverse = false;
				if(!cur->contig_nbs[left].ninfo.is_in())
				{
					has_reverse = true;
					cur->contig_nbs[left].ninfo.reverse();
				}
				order.push_back(vtuple(cur->contig_nbs[left].contigID, has_reverse,true));
				order.push_back(vtuple(cur->id, cur->contig_nbs[left].ninfo.right_isH(), false));
//             Don't have to reverse the right edge, even if it needs to be in logic
//				has_reverse = false;
//				if(cur->contig_nbs[right].ninfo.is_in()){
//					has_reverse = true;
//					cur->contig_nbs[right].ninfo.reverse();
//				}
				order.push_back(vtuple(cur->contig_nbs[right].contigID, cur->contig_nbs[right].ninfo.is_in(),true));
			}
			else
			{
				k_mer nb1 = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
				k_mer nb2  = cur->contig_nbs[0].nid;
				if(member_pos.find(nb1) != member_pos.end())
				{
					bool has_reverse = false;
					if(!cur->contig_nbs[0].ninfo.is_in())
					{
						has_reverse = true;
						cur->contig_nbs[0].ninfo.reverse();
					}
					order.push_back(vtuple(cur->contig_nbs[0].contigID, has_reverse,true));
					order.push_back(vtuple(cur->id, cur->contig_nbs[0].ninfo.right_isH(), false));
					next = nb1;
				}
				else
				{
					bool has_reverse = false;
					if(cur->contig_nbs[0].ninfo.is_in())
					{
						has_reverse = true;
						cur->contig_nbs[0].ninfo.reverse();
					}
					order.push_back(vtuple(cur->id, cur->contig_nbs[0].ninfo.left_isH(), false));
					order.push_back(vtuple(cur->contig_nbs[0].contigID, has_reverse,true));
					next = nb2;
				}
			}
		}
		//--- 2. process other vertices
		AmbiVertex * pre = cur;
		for(int i=1; i<count; i++)
		{
			cur = members[member_pos[next]];
			reordered[i] = cur;
			if(cur->type == V_1)
			{
				if(cur->ambi_nbs.size())
				{
					if(!cur->ambi_nbs[0].nb_info.is_in())
					{
						cur->ambi_nbs[0].nb_info.reverse();
					}
					order.push_back(vtuple(cur->id, cur->ambi_nbs[0].nb_info.right_isH(), false));
				}
				else
				{
					if(!cur->contig_nbs[0].ninfo.is_in())
					{
						cur->contig_nbs[0].ninfo.reverse();
					}
					order.push_back(vtuple(cur->id, cur->contig_nbs[0].ninfo.right_isH(), false));
				}
			}
			else
			{
				if(cur->ambi_nbs.size() == 2)
				{
					k_mer nb1 = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
					k_mer nb2 = get_neighbor(cur->id, cur->ambi_nbs[1].nb_info);
					int index;
					if(nb1 != pre->id)
					{
						index = 0;
						next = nb1;
					}
					else
					{
						index = 1;
						next = nb2;
					}
					if(cur->ambi_nbs[index].nb_info.is_in())
					{
						cur->ambi_nbs[index].nb_info.reverse();
					}
					order.push_back(vtuple(cur->id, cur->ambi_nbs[index].nb_info.left_isH(), false));
				}
				else if (cur->contig_nbs.size() == 2)
				{
					k_mer nb1 = cur->contig_nbs[0].nid;
					k_mer nb2 = cur->contig_nbs[1].nid;
					int index;
					if(nb1 != pre->id)
					{
						index = 0;
						next = nb1;
					}
					else
					{
						index = 1;
						next = nb2;
					}
					bool has_reverse = false;
					if(cur->contig_nbs[index].ninfo.is_in())
					{
						has_reverse = true;
						cur->contig_nbs[index].ninfo.reverse();
					}
					order.push_back(vtuple(cur->id, cur->contig_nbs[index].ninfo.left_isH(), false));
					order.push_back(vtuple(cur->contig_nbs[index].contigID,  has_reverse, true));
				}
				else
				{
					k_mer nb1 = get_neighbor(cur->id, cur->ambi_nbs[0].nb_info);
					k_mer nb2  = cur->contig_nbs[0].nid;
					if(nb1 != pre->id)
					{
						next = nb1;
						if(cur->ambi_nbs[0].nb_info.is_in())
						{
							cur->ambi_nbs[0].nb_info.reverse();
						}
						order.push_back(vtuple(cur->id, cur->ambi_nbs[0].nb_info.left_isH(), false));
					}
					else
					{
						next = nb2;
						bool has_reverse = false;
						if(cur->contig_nbs[0].ninfo.is_in())
						{
							has_reverse = true;
							cur->contig_nbs[0].ninfo.reverse();
						}
						order.push_back(vtuple(cur->id, cur->contig_nbs[0].ninfo.left_isH(), false));
						order.push_back(vtuple(cur->contig_nbs[0].contigID, has_reverse, true));
					}
				}
			}
			pre = cur;
		}
		members.swap(reordered);
		//end
	}

	int find(vector<MessageValue> & msgs, k_mer id)
	{
		for(int i = 0 ; i < msgs.size(); i++)
		{
			if(msgs[i].id == id)
				return i;
		}
		return -1;
	}

	void merge_sequence(AmbiContig * ambicontig, vector<MessageValue> & msgs, int index)
	{
		ATGC_bitmap tmp;
		if(order[index].iscontig)
		{
			int result = find(msgs, order[index].vid);
			if(result != -1)
			{
				if(order[index].isreverse)
				{
					tmp  = msgs[result].sequence.get_reverse();
				}
				else
				{
					tmp = msgs[result].sequence;
				}
			}
		}
		else
		{
			tmp.init(order[index].vid);
			if(order[index].isreverse)
			{
				tmp = tmp.get_reverse();
			}
		}
		if(index == 0)
		{
			for(int i = 0; i < tmp.length; i++)
			{
				ambicontig->seq.append(tmp.get(i));
			}
		}
		else
		{
			for(int i = mer_length - 1; i < tmp.length; i++)
			{
				ambicontig->seq.append(tmp.get(i));
			}
		}
	}

	AmbiContig * get_ambicontig(vector<MessageValue> & msgs)
	{
		AmbiContig * ambicontig = new AmbiContig;
		for(int i = 0; i < order.size(); i++)
		{
			merge_sequence(ambicontig, msgs, i);
		}
		return ambicontig;
	}
};

//AmbiVSet Key Order  =====================
//for defining hash_set<ContigVSet>
struct AmbiVSet_hash
{
	size_t operator()(AmbiVSet* x) const
	{
		return x->id;
	}
};

struct AmbiVSet_equal
{
public:
	bool operator()(AmbiVSet* a,  AmbiVSet* b) const
	{
		return a->id == b->id;
	}
};

class AmbiMergeWorker
{
public:
	typedef hash_set<AmbiVSet*, AmbiVSet_hash, AmbiVSet_equal> Groups;
	typedef Groups::iterator GroupIter;
	typedef vector<AmbiVSet*> GroupVector;
	typedef vector<AmbiContig*> AmbiContigVector;
	typedef vector<AmbiVertex *> AmbiVector;
	typedef vector<Contig *> ContigVector;

	int min_contig_length;
	DefaultHash<k_mer> hash;
	Groups groups;
	GroupVector groupVec;
	AmbiContigVector ambicontigs;
	AmbiVector ambiVec;
	ContigVector contigs;


	AmbiMergeWorker(int k, int min_contig_len)
	{
		set_mer_length(k);
		min_contig_length = min_contig_len;
	}

	~AmbiMergeWorker()
	{
		for(int i = 0; i < ambicontigs.size(); i++) delete ambicontigs[i];
	}

	void add_group(AmbiVSet* group)
	{
		GroupIter itr = groups.find(group);
		if(itr != groups.end())
		{
			AmbiVSet* cur = *itr;
			cur->merge(group);
		}
		else groups.insert(group);
	}

	void sync_groups()
	{
		//---------------------distribute groups to sending buffers----------------------------
		vector<GroupVector> _loaded_parts(_num_workers);
		for(GroupIter p = groups.begin(); p != groups.end(); p++)
		{
			AmbiVSet* cur = *p;
			_loaded_parts[hash(cur->id)].push_back(cur);
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
	}

	void groups_to_groupVec()
	{
		for(GroupIter p = groups.begin(); p != groups.end(); p++)  groupVec.push_back(*p);
		groups.clear();
	}

	void sync_contigs()
	{
		//---------------------distribute contigs to sending buffers----------------------------
		vector<ContigVector> _loaded_parts(_num_workers);
		for(ContigVector::iterator p = contigs.begin(); p != contigs.end(); p++)
		{
			Contig* cur = *p;
			_loaded_parts[hash(cur->id)].push_back(cur);
		}
		contigs.clear();
		//shuffle
		delete_after_all_to_all(_loaded_parts);
		//reset "contigs"
		for(int i = 0; i < _num_workers; i++)
		{
			ContigVector & vec = _loaded_parts[i];
			for(int j = 0; j < vec.size(); j++)
				contigs.push_back(vec[j]);
		}
	}

	void sync_AmbiVertex()
	{
		//---------------------distribute ambiVertex to sending buffers----------------------------
		vector<AmbiVector> _loaded_parts(_num_workers);
		for(AmbiVector::iterator p = ambiVec.begin(); p != ambiVec.end(); p++)
		{
			AmbiVertex* cur = *p;
			_loaded_parts[hash(cur->id)].push_back(cur);
		}
		ambiVec.clear();
		//shuffle
		delete_after_all_to_all(_loaded_parts);
		//reset "ambiVertex"
		for(int i = 0; i < _num_workers; i++)
		{
			AmbiVector & vec = _loaded_parts[i];
			for(int j = 0; j < vec.size(); j++)
				ambiVec.push_back(vec[j]);
		}
	}

	//==============================
	void add_vertex(char* line)
	{
		bool is_ambi = (line[strlen(line)-1] == '#');
		if(!is_ambi)
		{
			Contig * v = new Contig;
			v->parse(line);
			contigs.push_back(v);
		}
		else
		{
			AmbiVertex * v = new AmbiVertex;
			AmbiVSet* grp = new AmbiVSet;
			k_mer pred = v->parse(line);
			if(v->type == Vm_n)
				ambiVec.push_back(v);
			else
			{
				grp->id = pred;
				grp->add(v);
				add_group(grp);
			}
		}
	}

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

	void dump_partition(const char* outpath)
	{
		hdfsFS fs = getHdfsFS();
		BufferedWriter* writer = new BufferedWriter(outpath, fs, _my_rank);
		for(int i = 0; i < ambicontigs.size(); i++)
		{
			AmbiContig * cur = ambicontigs[i];
			if(cur->seq.length > min_contig_length)
				cur->dumpTo(writer, i);
		}
		delete writer;
		hdfsDisconnect(fs);
	}
	//=======================================================

	//Message passing between ambicontig, contig and ambiVertix
	void merge_ambicontig()
	{
		//   ------------------------------------Message Passing between ambiVSet and contigs---------------------------------------------------------------------
		DNAMessageBuffer<Contig, k_mer, k_mer, DefaultHash<k_mer> >* ambiVSet_message_buffer   = new DNAMessageBuffer<Contig, k_mer, k_mer, DefaultHash<k_mer> >;
		DNAMessageBuffer<AmbiVSet, k_mer, MessageValue, DefaultHash<k_mer> >* contigs_message_buffer  = new DNAMessageBuffer<AmbiVSet, k_mer, MessageValue, DefaultHash<k_mer> >;

		ambiVSet_message_buffer->init(contigs);
		contigs_message_buffer->init(groupVec);
		for(int i = 0; i < groupVec.size(); i++)
		{
			AmbiVSet* cur = groupVec[i];
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//add the following function in the case that we used SV-alg in previous step,
			//in that way, the group is not the head of this contig, we should find the head
#ifdef SV_USED
			k_mer head = cur->find_head_ambicontig();
			k_mer tail = cur->find_tail_ambicontig();
			if(head > tail) head = tail;
			cur->reorder_ambiVSet(head);
#else
			cur->reorder_ambiVSet();
#endif
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			for(int j = 0; j < cur->order.size(); j++)
			{
				if(cur->order[j].iscontig)
					ambiVSet_message_buffer->add_message(cur->order[j].vid, cur->id);
			}
		}
		ambiVSet_message_buffer->sync_messages();
		vector<vector<k_mer> >& v_msgbufs  = ambiVSet_message_buffer->get_v_msg_bufs();

		for(int i = 0 ; i < contigs.size(); i++)
		{
			if(v_msgbufs[i].size() != 0)
			{
				// every contig  will only be used once
				contigs[i]->mark = true;
				//only the first msg is needed
				//if there are more than one msg, some bugs must be triggered
				k_mer target = v_msgbufs[i][0];
				MessageValue msg;
				msg.id = contigs[i]->id;
				msg.sequence = contigs[i]->seq;
				contigs_message_buffer->add_message(target, msg);
				v_msgbufs[i].clear();
			}
			else
			{
				contigs[i]->mark = false;
			}
		}
		delete ambiVSet_message_buffer;

		contigs_message_buffer->sync_messages();
		vector<vector<MessageValue> >& con_msgbufs  = contigs_message_buffer->get_v_msg_bufs();
		for(int i = 0; i < groupVec.size(); i++)
		{
			AmbiVSet* cur = groupVec[i];
			AmbiContig * tmp = cur->get_ambicontig(con_msgbufs[i]);
			cur->free_members();
			delete cur;
			con_msgbufs[i].clear();
			ambicontigs.push_back(tmp);
		}
		delete contigs_message_buffer;
		//-------------------------------------- END--------------------------------------------------------

		//  ------------------------------Message Passing between ambiV and contigs---------------------------------
		DNAMessageBuffer<Contig, k_mer, k_mer, DefaultHash<k_mer> >* ambiV_message_buffer = new DNAMessageBuffer<Contig, k_mer, k_mer, DefaultHash<k_mer> >;
		DNAMessageBuffer<AmbiVertex, k_mer, MessageValue, DefaultHash<k_mer> >* remain_contigs_message_buffer = new DNAMessageBuffer<AmbiVertex, k_mer, MessageValue, DefaultHash<k_mer> >;

		ambiV_message_buffer->init(contigs);
		remain_contigs_message_buffer->init(ambiVec);
		for(int i = 0; i < ambiVec.size(); i++)
		{
			AmbiVertex* cur = ambiVec[i];
			for(int j = 0; j < cur->contig_nbs.size(); j++)
			{
				ambiV_message_buffer->add_message(cur->contig_nbs[j].contigID, cur->id);
			}
		}
		ambiV_message_buffer->sync_messages();
		vector<vector<k_mer> >& ambi_msgbufs  = ambiV_message_buffer->get_v_msg_bufs();
		for(int i = 0 ; i < contigs.size(); i++)
		{
			if(ambi_msgbufs[i].size() != 0 && contigs[i]->mark == false)
			{
				k_mer target = ambi_msgbufs[i][0];
				MessageValue msg;
				msg.id = contigs[i]->id;
				msg.sequence = contigs[i]->seq;
				remain_contigs_message_buffer->add_message(target, msg);
				ambi_msgbufs[i].clear();
			}
		}
		delete ambiV_message_buffer;

		remain_contigs_message_buffer->sync_messages();
		vector<vector<MessageValue> >& remain_msgbufs  = remain_contigs_message_buffer->get_v_msg_bufs();

		for(int i = 0; i < ambiVec.size(); i++)
		{
			for(int j = 0; j < remain_msgbufs[i].size(); j++)
			{
				AmbiContig * tmp = new AmbiContig;
				MessageValue & msg = remain_msgbufs[i][j];
				for(int k = 0; k < msg.sequence.length; k++)
					tmp->seq.append(msg.sequence.get(k));
				ambicontigs.push_back(tmp);
			}
			remain_msgbufs[i].clear();
			delete ambiVec[i];
		}
		delete remain_contigs_message_buffer;
		// --------------------------------------------------END----------------------------------------------
		// delete the loaded contigs, their information have been transmitted to ambicontigs
		for(int i = 0 ; i < contigs.size(); i++)
		{
			delete contigs[i];
		}
	}
	// ======================================================

	// run the worker
	void run(const MultiInputParams& params)
	{
		//check path + init
		if (_my_rank == MASTER_RANK)
		{
			if (dirCheck(params.input_paths, params.output_path.c_str(), _my_rank == MASTER_RANK, params.force_write) == -1)
				exit(-1);
		}
		init_timers();

		//dispatch splits
		ResetTimer(WORKER_TIMER);
		vector<vector<string> >* arrangement;
		if (_my_rank == MASTER_RANK)
		{
			arrangement = params.native_dispatcher ? dispatchLocality(params.input_paths) : dispatchRan(params.input_paths);
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
		sync_contigs();
		sync_AmbiVertex();

		StopTimer(WORKER_TIMER);
		PrintTimer("Sync Time", WORKER_TIMER);


		ResetTimer(WORKER_TIMER);
		// merge all the remained sequence into ambicontigs, deleted all other objects after message passing
		merge_ambicontig();
		//END

		StopTimer(WORKER_TIMER);
		PrintTimer("Merge Time", WORKER_TIMER);

		ResetTimer(WORKER_TIMER);

		dump_partition(params.output_path.c_str());

		StopTimer(WORKER_TIMER);
		PrintTimer("Dump Time", WORKER_TIMER);
	}
};

void AmbiMerge(string contig_path, string ambi_path, string out_path, int kmer, int min_contig_length)
{
	MultiInputParams ambmerge_p;
	ambmerge_p.input_paths.push_back(contig_path);
	ambmerge_p.input_paths.push_back(ambi_path);
	ambmerge_p.output_path = out_path;
	ambmerge_p.force_write=true;
	ambmerge_p.native_dispatcher=false;
	//	init_workers();
	AmbiMergeWorker ambmerge(kmer, min_contig_length);
	ambmerge.run(ambmerge_p);
	//	worker_finalize();
}

#endif /* AMBIMERGE_H_ */
