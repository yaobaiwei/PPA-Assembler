#ifndef SV_H_
#define SV_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"

using namespace std;

//if a k-mer is the end of a contig, append 1 before the k-mer
static bool v1_pol[4] = {false,false,true,true};
static bool v2_pol[4] = {false,true,false,true};

struct SVValue
{
	k_mer prev_D;
	k_mer D;
	vector<k_mer> neighbors;
	u8 type; //(1) 1, (2) 1-1, (3) m-n
	u32 bitmap;
	vector<u8> freqs;
};

ibinstream & operator<<(ibinstream & m, const SVValue & v)
{
	m<<v.prev_D;
	m<<v.D;
	m<<v.neighbors;
	m<<v.type;
	m<<v.bitmap;
	m<<v.freqs;
	return m;
}

obinstream & operator>>(obinstream & m, SVValue & v)
{
	m>>v.prev_D;
	m>>v.D;
	m>>v.neighbors;
	m>>v.type;
	m>>v.bitmap;
	m>>v.freqs;
	return m;
}

//====================================

int bitCount(u32 n)
{
	u32 tmp = n - ((n >>1) & 033333333333) - ((n >>2) & 011111111111);
	return ((tmp + (tmp >>3)) & 030707070707) % 63;
}

class SVVertex: public Vertex<k_mer, SVValue, k_mer>
{

public:
	inline k_mer is_contig_end(k_mer s) //return type treated as bool
	{
		return s >> (mer_length*2);
	}

	void get_neighbors(vector<k_mer> & collector)
	{
		int bound = value().type; //for type = 1 or 2
		if(value().type == 3) bound = bitCount(value().bitmap);
		//------
		for(int j = 0; j < 4; j ++)
		{
			int shift = getShift(v1_pol[j], v2_pol[j]);
			u32 shifted = (value().bitmap >> shift);
			//in-neighbors
			for(int i=0; i<8; i++)
			{
				if (shifted & ATGC_bits[i])
				{
					k_mer nb_id = get_neighbor(id, ATGC_bits[i], v1_pol[j], v2_pol[j]);
					collector.push_back(nb_id);
					if(collector.size() >= bound) return; //return earlier
				}
			}
		}
	}

	void get_neighbor_infos(vector<neighbor_info> & collector)
	{
		int bound = value().type; //for type = 1 or 2
		if(value().type == 3) bound = bitCount(value().bitmap);
		//------
		for(int j = 0; j < 4; j++)
		{
			int shift = getShift(v1_pol[j], v2_pol[j]);
			u32 shifted = (value().bitmap >> shift);
			//in-neighbors
			for(int i=0; i<8; i++)
			{
				if (shifted & ATGC_bits[i])
				{
					neighbor_info nb;
					nb.set_ATGC(bit2ATGC[i]);
					nb.set_in(i < 4);
					nb.set_left(v1_pol[j]);
					nb.set_right(v2_pol[j]);
					collector.push_back(nb);
					if(collector.size() >= bound) return; //return earlier
				}
			}
		}
	}

	void set_neighbors(MessageContainer & msgs) //assume ambiguous vertices have broadcasted msgs
	{
		//ambiguous vertices vote to halt; //unambiguous vertices set "preds" properly
		if(value().type == 1)
		{
			k_mer self = (id | END_MER);
			value().neighbors.push_back(self); //one pred is itself
			//another is the only neighbor
			if(msgs.size() == 1) value().neighbors.push_back(self);
			else
			{
				vector<k_mer> nbs;
				get_neighbors(nbs);
				value().neighbors.push_back(nbs[0]);
			}
		}
		else if(value().type == 2)
		{
			if(msgs.size() == 2) //itself is a contig
			{
				k_mer self = (id | END_MER);
				value().neighbors.push_back(self);
				value().neighbors.push_back(self);
			}
			else if(msgs.size() == 1)
			{
				k_mer self = (id | END_MER);
				value().neighbors.push_back(self);
				//find the other neighbor
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(msgs[0] == nbs[0]) value().neighbors.push_back(nbs[1]);
				else value().neighbors.push_back(nbs[0]);
			}
			else //msgs.size() == 0
			{
				get_neighbors(value().neighbors);
			}
		}
		else //value().type == 3
		{
			value().neighbors.swap(msgs); //collect ambi-nbs
		}
	}

	void treeInit_D()
	{
		for(int i = 0; i < value().neighbors.size(); i++)
		{
			k_mer nb = value().neighbors[i];
			if(!is_contig_end(nb) && nb < value().D)
			{
				value().D = nb;
			}
		}
	}

	void rtHook_1S()// = shortcut's request to w
	{
		// request to w
		k_mer Du=value().D;
		send_message(Du, id);
	}

	void rtHook_2R(MessageContainer & msgs)// = shortcut's respond by w
	{
		// respond by w
		k_mer Dw=value().D;
		for(int i=0; i<msgs.size(); i++)
		{
			k_mer requester=msgs[i];
			send_message(requester, Dw);
		}
	}

	void rtHook_2S()// = starhook's send D[v]
	{
		// send negated D[v]
		long long int Dv=value().D;
		for(int i = 0; i < value().neighbors.size(); i++)
		{
			k_mer nb = value().neighbors[i];
			if(!is_contig_end(nb))
			{
				send_message(nb, -Dv-1);//negate Dv to differentiate it from other msg types
			}
		}
	}

	void rtHook_3GDS(MessageContainer & msgs)//return whether a msg is sent
	{
		//set D[w]=min_v{D[v]} to allow fastest convergence, though any D[v] is ok (assuming (u, v) is accessed last)
		long long int Dw=-1;
		long long int Du=value().D;
		long long int Dv=-1;//pick the min
		for(int i=0; i<msgs.size(); i++)
		{
			long long int msg=msgs[i];
			if(msg>=0) Dw=msg;
			else// type==rtHook_2R_v
			{
				long long int cur=-msg-1;
				if(Dv==-1 || cur<Dv) Dv=cur;
			}
		}
		if(Dw==Du && Dv!=-1 && Dv<Du)//condition checking
		{
			send_message(Du, Dv);
		}
	}

	void rtHook_4GD(MessageContainer & msgs)// = starhook's write D[D[u]]
	{
		//set D[w]=min_v{D[v]} to allow fastest convergence, though any D[v] is ok (assuming (u, v) is accessed last)
		long long int Dv=-1;
		for(int i=0; i<msgs.size(); i++)
		{
			long long int cur=msgs[i];
			if(Dv==-1 || cur<Dv) Dv=cur;
		}
		if(Dv!=-1)
		{
			//value().pre_D = value().D;//YANDA: delete
			value().D =Dv;
		}
	}

	void shortcut_3GD(MessageContainer & msgs)  //D[u]=D[D[u]]
	{
		//value().pre_D = value().D;//YANDA: delete
		value().D =msgs[0];  //Once update the D[v], we should also update the Pre_D to keep the pre_D[v]
	}

	virtual void compute(MessageContainer & messages)
	{
		int cycle = 7;
		if(step_num() == 1)
		{
			if(value().type == 3)
			{
				vote_to_halt();
				vector<k_mer> nbs;
				get_neighbors(nbs);
				for(int i=0; i<nbs.size(); i++) send_message(nbs[i], id);
			}
		}
		else if(step_num() == 2)
		{
			set_neighbors(messages);
			if(value().type != 3)
			{
				treeInit_D();
				rtHook_1S();
			}
			else
			{
				vote_to_halt();
			}
		}
		else if(step_num() % cycle == 3)
		{
			if(value().type != 3)
			{
				rtHook_2R(messages);
				rtHook_2S();
			}
		}
		else if(step_num()% cycle == 4)
		{
			if(value().type != 3)
			{
				rtHook_3GDS(messages);
			}
		}
		else if(step_num()% cycle == 5)
		{
			if(value().type != 3)
			{
				rtHook_4GD(messages);
			}
		}
		else if(step_num()% cycle == 6)
		{
			if(value().type != 3)
			{
				rtHook_1S();
			}
		}
		else if(step_num()% cycle == 0)
		{
			if(value().type != 3)
			{
				rtHook_2R(messages);
			}
		}
		else if(step_num()% cycle == 1)
		{
			if(value().type != 3)
			{
				shortcut_3GD(messages);
			}
		}
		else if(step_num() % cycle == 2)
		{
			bool* agg=(bool*)getAgg();
			if(*agg)
			{
				vote_to_halt();
				return;
			}
			rtHook_1S();
		}
	}
};

//====================================
u8 get_type(u32 bitmap)
{
	int num = bitCount(bitmap);
	if(num == 1) return 1;
	else if(num == 2)
	{
		int count[2] = {0,0}; // inH/outL, outH/inL
		//LL: inL, outL
		//LH: inH, outL
		//HL: inL, outH
		//HH: inH, outH
		int pos[8] = {0,1,0,0,1,1,1,0};
		int cur = 0;
		for(int j = 0; j < 4; j++)
		{
			for(int i = 0; i < 4; i++)
			{
				if((bitmap >> i) & 1u)
				{
					int & cnt = count[pos[cur]];
					cnt ++;
					if(cnt > 1) return 3;
				}
			}
			cur++;
			//----
			for(int i = 4; i < 8; i++)
			{
				if((bitmap >> i) & 1u)
				{
					int & cnt = count[pos[cur]];
					cnt ++;
					if(cnt > 1) return 3;
				}
			}
			cur++;
			//----
			bitmap >>= 8;
		}
		return 2;
	}
	else return 3; //0 or >2
}

class SVAgg:public Aggregator<SVVertex, bool, bool>
{
private:
	bool AND;

public:
	SVAgg()
	{
		AND = true;
	}

	virtual void init()
	{
		AND=true;
	}

	virtual void stepPartial(SVVertex* v)
	{
		if(step_num() % 7 == 1 && step_num() > 1)
			if(v->value().prev_D != v->value().D)
			{
				AND=false;
				v->value().prev_D = v->value().D;
			}
	}

	virtual void stepFinal(bool* part)
	{
		if(*part==false) AND=false;
	}

	virtual bool* finishPartial()
	{
		return &AND;
	}
	virtual bool* finishFinal()
	{
		return &AND;
	}
};

class SVWorker:public Worker<SVVertex, SVAgg>
{
	char buf[100];

public:

	SVWorker(int k)
	{
		set_mer_length(k);
		END_MER = 1ull << (2*mer_length);
	}

	virtual SVVertex* toVertex(char* line)
	{
		char * pch;
		SVVertex* v=new SVVertex;
		pch=strtok(line, "\t");
		v->id = strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		u32 bitmap = (u32)strtoul(pch, NULL, 10);
		v->value().bitmap = bitmap;
		pch=strtok(NULL, " ");
		int size = atoi(pch);
		for(int i = 0; i < size; i++)
		{
			pch=strtok(NULL, " ");
			v->value().freqs.push_back((u8)atoi(pch));
		}
		v->value().type = get_type(bitmap);
		v->value().prev_D = v->id;
		v->value().D = v->id;
		return v;
	}

	virtual void toline(SVVertex* v, vector<BufferedWriter *> & writers)
	{
		//if type = 3, output to amb_out: vid \t num_nbs nb1 freq1 (nb2 freq2) ... //only ambi-neighbors
		//otherwise, line format: vid \t smaller_pred type(i.e. num_nbs) nb1 freq1 (nb2 freq2) ...
		SVValue & val = v->value();
		if(val.type == 3)
		{
			sprintf(buf, "%llu\t", v->id);
//			convert(v->id, buf, mer_length);  // human read
//			buf[mer_length] = '\0';
			writers[1]->write(buf);
			//------
			//build ID -> pos map
			vector<k_mer> nbs;
			v->get_neighbors(nbs);
			hash_map<k_mer, int> nb_pos;
			for(int i = 0; i < nbs.size(); i++ )
			{
				nb_pos[nbs[i]] = i;
			}
			//------
			vector<neighbor_info> nb_infos;
			v->get_neighbor_infos(nb_infos);
			vector<u32> counts;
			parse_vints(counts, val.freqs);
			//------
			sprintf(buf, "%u",val.neighbors.size());
			writers[1]->write(buf);
			for(int i = 0; i < val.neighbors.size(); i++)
			{
				int pos = nb_pos[val.neighbors[i]]; //get pos of ambi-nb
				sprintf(buf, " %u %u", nb_infos[pos].bitmap, counts[pos]);
				writers[1]->write(buf);
			}
			writers[1]->write(" #\n"); //# is a special tag indicating that it's an ambi-node
		}
		else
		{
			vector<neighbor_info> nbs;
			v->get_neighbor_infos(nbs);
			//----
			vector<u32> counts;
			parse_vints(counts, val.freqs);
			//----
			sprintf(buf, "%llu\t%llu %u", v->id, v->value().D, counts.size());
			writers[0]->write(buf);
			for(int i=0; i<counts.size(); i++)
			{
				sprintf(buf, " %u %u", nbs[i].bitmap, counts[i]);
				writers[0]->write(buf);
			}
			writers[0]->write("\n");
		}
	}
};

void SV(string in_path, string unamb_path, string amb_path, int length)
{
	MultiOutputParams param;
	param.input_path=in_path;
	param.output_paths.push_back(unamb_path);
	param.output_paths.push_back(amb_path);
	param.force_write=true;
	param.native_dispatcher=false;
	//	init_workers();
	SVWorker worker(length);
	SVAgg agg = SVAgg();
	worker.setAggregator(&agg);
	worker.run(param);
	//	worker_finalize();
}


#endif /* SV_H_ */
