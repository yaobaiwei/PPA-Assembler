#ifndef LISTRANK_H
#define LISTRANK_H

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"

using namespace std;

//if a k-mer is the end of a contig, append 1 before the k-mer
static bool v1_pol[4] = {false,false,true,true};
static bool v2_pol[4] = {false,true,false,true};
bool is_SV = false;

struct LRValue
{
	vector<k_mer> preds; //two preds in two directions //or for ambi-node, the ambi-nbs
	u8 type; //(1) 1, (2) 1-1, (3) m-n
	u32 bitmap;
	vector<u8> freqs;
};

ibinstream & operator<<(ibinstream & m, const LRValue & v)
{
	m<<v.preds;
	m<<v.type;
	m<<v.bitmap;
	m<<v.freqs;
	return m;
}

obinstream & operator>>(obinstream & m, LRValue & v)
{
	m>>v.preds;
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

class LRVertex: public Vertex<k_mer, LRValue, k_mer>
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

	void set_preds(MessageContainer & msgs) //assume ambiguous vertices have broadcasted msgs
	{
		//ambiguous vertices vote to halt; //unambiguous vertices set "preds" properly
		if(value().type == 1)
		{
			k_mer self = (id | END_MER);
			value().preds.push_back(self); //one pred is itself
			//another is the only neighbor
			if(msgs.size() == 1) value().preds.push_back(self);
			else
			{
				vector<k_mer> nbs;
				get_neighbors(nbs);
				value().preds.push_back(nbs[0]);
			}
		}
		else if(value().type == 2)
		{
			if(msgs.size() == 2) //itself is a contig
			{
				k_mer self = (id | END_MER);
				value().preds.push_back(self);
				value().preds.push_back(self);
			}
			else if(msgs.size() == 1)
			{
				k_mer self = (id | END_MER);
				value().preds.push_back(self);
				//find the other neighbor
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(msgs[0] == nbs[0]) value().preds.push_back(nbs[1]);
				else value().preds.push_back(nbs[0]);
			}
			else //msgs.size() == 0
			{
				get_neighbors(value().preds);
			}
		}
		else //value().type == 3
		{
			value().preds.swap(msgs); //collect ambi-nbs
		}
	}

	void LR_req()
	{
		//send myself as a requester to both preds
		if(!is_contig_end(value().preds[0])) send_message(value().preds[0], id);
		if(!is_contig_end(value().preds[1])) send_message(value().preds[1], id);
	}

	void LR_resp(MessageContainer & msgs)
	{
		for(int i=0; i<msgs.size(); i++)
		{
			if(msgs[i] == (value().preds[0] & kick)) send_message(msgs[i], value().preds[1]);
			else send_message(msgs[i], value().preds[0]);
		}
	}

	void treeInit_D()
	{
		//set D[u]=min{v} to allow fastest convergence, though any v is ok (assuming (u, v) is accessed last)
		vector<k_mer> nbs;
		get_neighbors(nbs);
		for(int i=0; i<nbs.size(); i++)
		{
			k_mer nb= nbs[i];
			if(nb<value().preds[1]) value().preds[1]=nb;
		}
	}

	void rtHook_1S()// = shortcut's request to w
	{
		// request to w
		k_mer Du=value().preds[1];
		send_message(Du, id);
	}

	void rtHook_2R(MessageContainer & msgs)// = shortcut's respond by w
	{
		// respond by w
		k_mer Dw=value().preds[1];
		for(int i=0; i<msgs.size(); i++)
		{
			k_mer requester=msgs[i];
			send_message(requester, Dw);
		}
	}

	void rtHook_2S()// = starhook's send D[v]
	{
		// send negated D[v]
		long long int Dv=value().preds[1];
		vector<k_mer> nbs;
		get_neighbors(nbs);
		for(int i=0; i<nbs.size(); i++)
		{
			k_mer nb=nbs[i];
			send_message(nb, -Dv-1);//negate Dv to differentiate it from other msg types
		}
	}//in fact, a combiner with MIN operator can be used here

	void rtHook_3GDS(MessageContainer & msgs)//return whether a msg is sent
	{
		//set D[w]=min_v{D[v]} to allow fastest convergence, though any D[v] is ok (assuming (u, v) is accessed last)
		long long int Dw=-1;
		long long int Du=value().preds[1];
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
			value().preds[1] =Dv;
		}
	}

	void shortcut_3GD(MessageContainer & msgs)  //D[u]=D[D[u]]
	{
		//value().pre_D = value().D;//YANDA: delete
		value().preds[1] =msgs[0];  //Once update the D[v], we should also update the Pre_D to keep the pre_D[v]
	}

	virtual void compute(MessageContainer & messages)
	{
		if(! is_SV)
		{
			if(step_num() == 1)
			{
				//ambiguous vertices send messages to neighbors
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
				vote_to_halt();
				set_preds(messages);
				if(value().type != 3) LR_req(); //otherwise, there is no preds
			}
			else if(step_num() % 2 == 1)
			{
				vote_to_halt();
				LR_resp(messages);
			}
			else if(step_num() % 2 == 0)
			{
				vote_to_halt();
				//update preds
				if(is_contig_end(value().preds[0])) messages.push_back(value().preds[0]);
				if(is_contig_end(value().preds[1])) messages.push_back(value().preds[1]);
				value().preds.swap(messages);
				//---
				LR_req();
			}
		}
		else
		{
			int cycle = 7;
			if(step_num() == 1)
			{
				value().preds[0] = value().preds[1] = id;
				treeInit_D();
				rtHook_1S();
			}
			else if(step_num() % cycle == 2)
			{
				rtHook_2R(messages);
				rtHook_2S();
			}
			else if(step_num() % cycle == 3)
			{
				rtHook_3GDS(messages);
			}
			else if(step_num() % cycle == 4)
			{
				rtHook_4GD(messages);
			}

			else if(step_num() % cycle == 5)
			{
				rtHook_1S();
			}
			else if(step_num() % cycle == 6)
			{
				rtHook_2R(messages);
			}
			else if(step_num() % cycle == 0)
			{
				shortcut_3GD(messages);
			}

			else if(step_num() % cycle == 1)
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

class LRAgg:public Aggregator<LRVertex, bool, bool>
{
private:
	bool AND;
	long long int msgs_size;

public:
	LRAgg(long long int msgs)
	{
		msgs_size = msgs;
		AND = true;
	}

	virtual void init()
	{
		AND=true;
	}

	virtual void stepPartial(LRVertex* v)
	{
		if(is_SV)
		{
			if (step_num() % 7 == 0)
				if (v->value().preds[0]  != v->value().preds[1])
				{
					AND = false;
					v->value().preds[0] = v->value().preds[1];
				}
		}
	}

	virtual void stepFinal(bool* part)
	{
		if(*part==false) AND=false;
	}

	virtual bool* finishPartial()
	{
		if(!is_SV)
		{
			if(step_num() % 2 == 1)
			{
				if(msgs_size != get_step_msg_num())
					msgs_size = get_step_msg_num();
				else
				{
					is_SV = true;
					global_step_num = 0;
				}
			}
		}
		return &AND;
	}

	virtual bool* finishFinal()
	{
		if(!is_SV)
		{
			if(step_num() % 2 == 1)
			{
				if(msgs_size != get_step_msg_num())
					msgs_size = get_step_msg_num();
				else
				{
					is_SV = true;
					global_step_num = 0;
				}
			}
		}
		return &AND;
	}
};

class LRWorker:public Worker<LRVertex, LRAgg>
{
	char buf[100];

public:

	LRWorker(int k)
	{
		set_mer_length(k);
		END_MER = 1ull << (2*mer_length);
	}

	virtual LRVertex* toVertex(char* line)
	{
		char * pch;
		LRVertex* v=new LRVertex;
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
		return v;
	}

	virtual void toline(LRVertex* v, vector<BufferedWriter *> & writers)
	{
		//if type = 3, output to amb_out: vid \t num_nbs nb1 freq1 (nb2 freq2) ... //only ambi-neighbors
		//otherwise, line format: vid \t smaller_pred type(i.e. num_nbs) nb1 freq1 (nb2 freq2) ...
		LRValue & val = v->value();
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
			sprintf(buf, "%u",val.preds.size());
			writers[1]->write(buf);
			for(int i = 0; i < val.preds.size(); i++)
			{
				int pos = nb_pos[val.preds[i]]; //get pos of ambi-nb
				sprintf(buf, " %u %u", nb_infos[pos].bitmap, counts[pos]);
				writers[1]->write(buf);
			}
			writers[1]->write(" #\n"); //# is a special tag indicating that it's an ambi-node
		}
		else
		{
			k_mer pred1 = (val.preds[0] & kick);
			k_mer pred2 = (val.preds[1] & kick);
			if(pred1 > pred2)
			{
				pred1 = pred2;
			}
			//----
			vector<neighbor_info> nbs;
			v->get_neighbor_infos(nbs);
			//----
			vector<u32> counts;
			parse_vints(counts, val.freqs);
			//----
			sprintf(buf, "%llu\t%llu %u", v->id, pred1, counts.size());
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

void ListRank(string in_path, string unamb_path, string amb_path, int length)
{
	MultiOutputParams param;
	param.input_path=in_path;
	param.output_paths.push_back(unamb_path);
	param.output_paths.push_back(amb_path);
	param.force_write=true;
	param.native_dispatcher=false;
	//	init_workers();
	LRWorker worker(length);
	LRAgg agg = LRAgg(-1);
	worker.setAggregator(&agg);
	worker.run(param);
	//	worker_finalize();
}

#endif
