#ifndef AMBICONNECT_H_
#define AMBICONNECT_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"
using namespace std;

struct ConnContigValue
{
	k_mer in_neighbor;
	bool in_pol;
	u32 in_count;
	k_mer out_neighbor;
	bool out_pol;
	u32 out_count;
	int freq;
	ATGC_bitmap seq;
};

ibinstream & operator<<(ibinstream & m, const ConnContigValue & v)
{
	m << v.in_neighbor;
	m << v.in_pol;
	m << v.in_count;
	m << v.out_neighbor;
	m << v.out_pol;
	m << v.out_count;
	m << v.freq;
	m << v.seq;
	return m;
}

obinstream & operator>>(obinstream & m,  ConnContigValue & v)
{
	m >> v.in_neighbor;
	m >> v.in_pol;
	m >> v.in_count;
	m >> v.out_neighbor;
	m >> v.out_pol;
	m >> v.out_count;
	m >> v.freq;
	m >> v.seq;
	return m;
}

struct ConnAmbiValue
{
	vector<neighbor_info> ambi_nbs;
	vector<u32> freqs;
	vector<ContigNB> contig_nbs;
};

ibinstream & operator<<(ibinstream & m, const ConnAmbiValue & v)
{
	m << v.ambi_nbs;
	m << v.freqs;
	m << v.contig_nbs;
	return m;
}

obinstream & operator>>(obinstream & m,  ConnAmbiValue & v)
{
	m >> v.ambi_nbs;
	m >> v.freqs;
	m >> v.contig_nbs;
	return m;
}

//=============================
static const u8 AMBI_TYPE = 0;
static const u8 CONTIG_TYPE = 1; //non-tip contig
static const u8 TIP_TYPE = 2; //tip contig
//=============================

struct ConnValue
{
	u8 type;
	void * value;

	~ConnValue()
	{
		if(type == AMBI_TYPE)
		{
			ConnAmbiValue * val = (ConnAmbiValue *)value;
			delete val;
		}
		else
		{
			ConnContigValue * val = (ConnContigValue *)value;
			delete val;
		}
	}
};

ibinstream & operator<<(ibinstream & m, const ConnValue & v)
{
	m << v.type;
	if(v.type == AMBI_TYPE) m << *((ConnAmbiValue*)v.value);
	else m << *((ConnContigValue*)v.value);
	return m;
}

obinstream & operator>>(obinstream & m, ConnValue & v)
{
	m >> v.type;
	if(v.type == AMBI_TYPE)
	{
		v.value = new ConnAmbiValue;
		m >> *((ConnAmbiValue*)v.value);
	}
	else
	{
		v.value = new ConnContigValue;
		m >> *((ConnContigValue*)v.value);
	}
	return m;
}

static int tipLength_threshold;

class ConnVertex: public Vertex<k_mer, ConnValue, ContigNB>
{
public:
	u8 getType(ConnContigValue * value)
	{
		if(value->seq.length > tipLength_threshold) return CONTIG_TYPE;
		if(value->in_neighbor == NULL_MER) return TIP_TYPE;
		if(value->out_neighbor == NULL_MER) return TIP_TYPE;
		return CONTIG_TYPE;
	}

	virtual void compute(MessageContainer & messages)
	{
		if(step_num() == 1)
		{
			if(value().type == CONTIG_TYPE)
			{
				ConnContigValue * val = (ConnContigValue*)(value().value);
				int type = getType(val);
				value().type = type;
				if(type == CONTIG_TYPE)
				{
					ContigNB msg;
					msg.contigID = id;
					msg.length = val->seq.length;
					msg.ninfo.set_left(val->in_pol);
					msg.ninfo.set_right(val->out_pol);
					if(val->in_neighbor != NULL_MER)
					{
						msg.ninfo.set_in(false);
						msg.nid = val->out_neighbor;
						send_message(val->in_neighbor, msg);
					}
					if(val->out_neighbor != NULL_MER)
					{
						msg.ninfo.set_in(true);
						msg.nid = val->in_neighbor;
						send_message(val->out_neighbor, msg);
					}
				}
			}
		}
		else if(step_num() == 2)
		{
			if(value().type == AMBI_TYPE)
			{
				ConnAmbiValue * val =  (ConnAmbiValue*)(value().value);
				int size = messages.size();
				for(int i = 0; i < size; i++)
				{
					val->contig_nbs.push_back(messages[i]);
				}
			}
		}
		vote_to_halt();
	}
};

class ConnectWorker:public Worker<ConnVertex>
{
	char buf[1024];

public:
	ConnectWorker(int k, int tip_len)
	{
		set_mer_length(k);
		tipLength_threshold = k - 1 + tip_len;
	}

	virtual ConnVertex* toVertex(char* line)
	{
		bool is_ambi = (line[strlen(line)-1] == '#');
		char * pch;
		ConnVertex* v=new ConnVertex;
		pch=strtok(line, "\t");
		v->id = strtoull(pch, NULL, 10);
		if(is_ambi)
		{
			v->value().type = AMBI_TYPE;
			ConnAmbiValue * val = new ConnAmbiValue;
			pch=strtok(NULL, " ");
			int size = atoi(pch);
			for(int i = 0; i< size; i++)
			{
				pch=strtok(NULL, " ");
				neighbor_info nb;
				nb.bitmap = (u8)atoi(pch);
				val->ambi_nbs.push_back(nb);
				pch=strtok(NULL, " ");
				val->freqs.push_back(atoi(pch));
			}
			v->value().value = val;
		}
		else
		{
			v->value().type = CONTIG_TYPE;
			ConnContigValue * val = new ConnContigValue;
			pch=strtok(NULL, " ");
			val->in_neighbor= strtoull(pch, NULL, 10);
			pch=strtok(NULL, " ");
			val->in_pol = atoi(pch);
			pch=strtok(NULL, " ");
			val->in_count = atoi(pch);
			pch=strtok(NULL, " ");
			val->out_neighbor = strtoull(pch, NULL, 10);
			pch=strtok(NULL, " ");
			val->out_pol = atoi(pch);
			pch=strtok(NULL, " ");
			val->out_count = atoi(pch);
			pch=strtok(NULL, " ");
			val->freq = atoi(pch);
			pch=strtok(NULL, " ");
			val->seq.length = atoi(pch);
			pch=strtok(NULL, " ");
			int size  =  atoi(pch);
			for(int i = 0; i < size; i++)
			{
				pch=strtok(NULL, " ");
				val->seq.seq.push_back((u8)atoi(pch));
			}
			v->value().value = val;
		}
		return v;
	}

	virtual void toline(ConnVertex* v, BufferedWriter& writer)
	{
		ConnValue & value = v->value();
		if(value.type == AMBI_TYPE)
		{
			ConnAmbiValue * val = (ConnAmbiValue *)(value.value);
			int size_amb = val->ambi_nbs.size();
			int size_con = val->contig_nbs.size();
			if(size_amb > 0 || size_con > 0)
			{
				sprintf(buf, "%llu\t%d", v->id, size_amb);
				writer.write(buf);
				for(int i=0; i<size_amb; i++)
				{
					sprintf(buf, " %u %u", val->ambi_nbs[i].bitmap, val->freqs[i]);
					writer.write(buf);
				}
				sprintf(buf, " %d", size_con);
				writer.write(buf);
				for(int i=0; i<size_con; i++)
				{
					ContigNB & cur = val->contig_nbs[i];
					sprintf(buf, " %llu %u %llu %d", cur.nid, cur.ninfo.bitmap, cur.contigID, cur.length);
					writer.write(buf);
				}
				writer.write("\n");
			}
		}
	}
};

void AmbiConnect(string amb_path, string contig_path, string out_path, int mer_len, int tip_len)
{
	MultiInputParams param;
	param.input_paths.push_back(amb_path);
	param.input_paths.push_back(contig_path);
	param.output_path = out_path;
	param.force_write=true;
	param.native_dispatcher=false;
//	init_workers();
	ConnectWorker worker(mer_len, tip_len);
	worker.run(param);
//	worker_finalize();
}

#endif
