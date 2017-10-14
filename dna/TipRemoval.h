#ifndef TIPREMOVAL_H_
#define TIPREMOVAL_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"

using namespace std;

static int tipLen_threshold;

//=================================

//v-status:
static const u8 Deleted = 4;
static const u8 Finished = 5;
static const u8 Normal = 6;

struct TRVertexValue
{
	u8 status;
	u8 type;
	vector<AmbiNB> ambi_nbs;
	vector<ContigNB> contig_nbs;
};

ibinstream & operator<<(ibinstream & m, const TRVertexValue & v)
{
	m << v.status;
	m << v.type;
	m << v.ambi_nbs;
	m << v.contig_nbs;
	return m;
}

obinstream & operator>>(obinstream & m,  TRVertexValue & v)
{
	m >> v.status;
	m >> v.type;
	m >> v.ambi_nbs;
	m >> v.contig_nbs;
	return m;
}

//=================================

struct TRMsg
{
	k_mer source;
	int length;
	char status;
};

ibinstream & operator<<(ibinstream & m, const TRMsg & v)
{
	m << v.source;
	m << v.length;
	m << v.status;
	return m;
}

obinstream & operator>>(obinstream & m,  TRMsg & v)
{
	m >> v.source;
	m >> v.length;
	m >> v.status;
	return m;
}

//=================================

class TRVertex: public Vertex<k_mer, TRVertexValue, TRMsg>
{
public:
	int getType()
	{
		int size = value().ambi_nbs.size() + value().contig_nbs.size();
		if(size == 1) return V_1;
		else if(size == 2)
		{
			int count[2] = {0}; // inH/outL, outH/inL
			for(int j = 0; j < value().ambi_nbs.size(); j++)
			{
				AmbiNB & nb = value().ambi_nbs[j];
				neighbor_info & ninfo = nb.nb_info;
				if((ninfo.is_in() && ninfo.right_isH()) || (!ninfo.is_in() && !ninfo.left_isH()))
				{
					count[0]++;
					if(count[0] > 1) return Vm_n;
				}
				else
				{
					count[1]++;
					if(count[1] > 1) return Vm_n;
				}
			}
			for(int j = 0 ; j < value().contig_nbs.size(); j++)
			{
				ContigNB & nb = value().contig_nbs[j];
				neighbor_info & ninfo = nb.ninfo;
				if((ninfo.is_in() && ninfo.right_isH()) || (!ninfo.is_in() && !ninfo.left_isH()))
				{
					count[0]++;
					if(count[0] > 1) return Vm_n;
				}
				else
				{
					count[1]++;
					if(count[1] > 1) return Vm_n;
				}
			}
			return V1_1;
		}
		else return Vm_n; //0 or >2
	}

	int forward_length(int msg_length = 0, k_mer msg_src = NULL_MER)
	{
		for(int i = 0; i < value().ambi_nbs.size(); i++)
		{
			AmbiNB & cur = value().ambi_nbs[i];
			k_mer target = get_neighbor(id, cur.nb_info);
			if(msg_src != target)
			{
				TRMsg msg;
				msg.source = id;
				msg.length = msg_length + 1;
				msg.status = Normal;
				send_message(target, msg);
				return -1;
			}
		}
		//------
		for(int i = 0; i < value().contig_nbs.size(); i++)
		{
			ContigNB & cur = value().contig_nbs[i];
			if(msg_src != cur.nid)
			{
				if(cur.nid == NULL_MER) return msg_length + 1 +cur.length;
				TRMsg msg;
				msg.source = id;
				msg.length = msg_length + 1 + cur.length;
				msg.status = Normal;
				send_message(cur.nid, msg);
				return -1;
			}
		}
		return value().contig_nbs[0].length + 1;
	}

	void forward_status(u8 status, k_mer msg_src)
	{
		for(int i = 0; i < value().ambi_nbs.size(); i++)
		{
			AmbiNB & cur = value().ambi_nbs[i];
			k_mer target = get_neighbor(id, cur.nb_info);
			if( msg_src != target)
			{
				TRMsg msg;
				msg.source = id;
				msg.length = 0;
				msg.status = status;
				send_message(target, msg);
				return;
			}
		}
		//------
		for(int i = 0; i < value().contig_nbs.size(); i++)
		{
			ContigNB & cur = value().contig_nbs[i];
			if(msg_src != cur.nid)
			{
				TRMsg msg;
				msg.source = id;
				msg.length = 0;
				msg.status = status;
				send_message(cur.nid, msg);
				return;
			}
		}
	}

	void respond(k_mer target, u8 status, int len = 0)
	{
		TRMsg msg;
		msg.source = id;
		msg.length = len;
		msg.status = status;
		send_message(target, msg);
	}

	void delete_edge(k_mer nb)
	{
		vector<AmbiNB>::iterator it;
		for(it = value().ambi_nbs.begin(); it != value().ambi_nbs.end(); it++)
		{
			k_mer target = get_neighbor(id, it->nb_info);
			if(target == nb)
			{
				value().ambi_nbs.erase(it);
				return;
			}
		}
		vector<ContigNB>::iterator itr;
		for(itr = value().contig_nbs.begin(); itr != value().contig_nbs.end(); itr++)
		{
			if(itr->nid == nb)
			{
				value().contig_nbs.erase(itr);
				return;
			}
		}
	}

	virtual void compute(MessageContainer & messages)
	{
		vote_to_halt();
		if(step_num() == 1)
		{
			if(value().type == V_1 && value().status == Normal)
			{
				int result = forward_length();
				if(result != -1)
				{
					if(result <= tipLen_threshold)
						value().status = Deleted;
					else
						value().status = Finished;
				}
			}
		}
		else
		{
			if(value().status != Normal)
			{
				for(int i = 0; i < messages.size(); i++)
				{
					if(messages[i].length != -1)
						respond(messages[i].source, value().status, -1);
				}
				return;
			}
			//------
			int type = value().type;
			if(type == V_1)
			{
				if(messages[0].status != Normal || messages[0].length == -1)
				{
					value().status = messages[0].status;
					return;
				}
				//V1 - {V1-1} - V1
				if(messages[0].length + 1 <= tipLen_threshold)
				{
					value().status = Deleted;
					respond(messages[0].source, Deleted);
				}
				else
				{
					value().status = Finished;
					respond(messages[0].source, Finished);
				}
			}
			else if(type == V1_1)
			{
				for(int i = 0; i < messages.size(); i++)
				{
					if(messages[i].status != Normal || messages[i].length == -1)
					{
						value().status = messages[i].status;
						forward_status(messages[i].status, messages[i].source);
					}
					else
					{
						int result = forward_length(messages[i].length, messages[i].source);
						if(result != -1)
						{
							if(result <= tipLen_threshold)
							{
								value().status = Deleted;
								forward_status(Deleted, NULL_MER);
							}
							else
							{
								value().status = Finished;
								forward_status(Finished, NULL_MER);
							}
						}
					}
				}
			}
			else if(type == Vm_n)
			{
				for(int i = 0; i < messages.size(); i++)
				{
					if(messages[i].length + 1 <= tipLen_threshold)
					{
						respond(messages[i].source, Deleted);
						delete_edge(messages[i].source);
					}
					else respond(messages[i].source, Finished);
				}
			}
		}
	}

};

class TRWorker:public Worker<TRVertex>
{
	char buf[1024];

public:

	TRWorker(int k, int tip_len)
	{
		set_mer_length(k);
		tipLen_threshold = k - 1 + tip_len;
	}

	virtual TRVertex* toVertex(char* line)
	{
		char * pch;
		TRVertex* v=new TRVertex;
		pch=strtok(line, "\t");
		v->id = strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		int size = atoi(pch);
		for(int i = 0; i < size; i++)
		{
			pch=strtok(NULL, " ");
			AmbiNB ambiNB;
			ambiNB.nb_info.bitmap = (u8)atoi(pch);
			pch=strtok(NULL, " ");
			ambiNB.count = atoi(pch);
			v->value().ambi_nbs.push_back(ambiNB);
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
			v->value().contig_nbs.push_back(contigNB);
		}
		v->value().type = v->getType();
		v->value().status = Normal;
		return v;
	}

	virtual void toline(TRVertex* v, BufferedWriter& writer)
	{
		TRVertexValue & value = v->value();
		if(value.status != Deleted)
		{
			int size = value.ambi_nbs.size();
			sprintf(buf, "%llu\t%u %d", v->id, value.type, size);
			writer.write(buf);
			//------
			for(int i = 0; i < size; i++)
			{
				AmbiNB & nb = value.ambi_nbs[i];
				sprintf(buf, " %u %d", nb.nb_info.bitmap, nb.count);
				writer.write(buf);
			}
			//------
			size = value.contig_nbs.size();
			sprintf(buf, " %d", size);
			writer.write(buf);
			for(int i = 0; i < size; i++)
			{
				ContigNB & nb = value.contig_nbs[i];
				sprintf(buf, " %llu %u %llu %d", nb.nid, nb.ninfo.bitmap, nb.contigID, nb.length);
				writer.write(buf);
			}
			writer.write("\n");
		}
	}

	virtual bool phase_continue(vector<TRVertex*> & vertexes) //this is what user specifies!!!!!!
	{
		if(phase_num() == 1) return true;
		else
		{
			bool has_v1 = false;
			for(int i = 0; i <  vertexes.size(); i++)
			{
				TRVertex * cur = vertexes[i];
				if(cur->value().status == Normal)
				{
					int type = cur->getType();
					cur->value().type = type;
					if(type == V_1) has_v1 = true;
				}
			}
			return all_bor(has_v1);
		}
	}
};

void TipRemoval(string amb_path, string out_path, int mer_len, int tip_len)
{
	WorkerParams param;
	param.input_path = amb_path;
	param.output_path = out_path;
	param.force_write=true;
	param.native_dispatcher=false;
//	init_workers();
	TRWorker worker(mer_len, tip_len);
	worker.runPhase(param);
//	worker_finalize();
}

#endif
