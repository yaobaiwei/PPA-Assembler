#ifndef AMBILISTRANK_H_
#define AMBILISTRANK_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"

using namespace std;

bool Amb_is_SV = false;

struct AmbLRValue
{
	u8 type;
	vector<k_mer> preds;
	vector<AmbiNB> ambi_nbs;
	vector<ContigNB> contig_nbs;
};

ibinstream & operator<<(ibinstream & m, const AmbLRValue & v)
{
	m << v.type;
	m << v.preds;
	m << v.ambi_nbs;
	m << v.contig_nbs;
	return m;
}

obinstream & operator>>(obinstream & m,  AmbLRValue & v)
{
	m >> v.type;
	m >> v.preds;
	m >> v.ambi_nbs;
	m >> v.contig_nbs;
	return m;
}

class AmbLRVertex: public Vertex<k_mer, AmbLRValue, k_mer>
{
public:

	inline k_mer is_contig_end(k_mer s) //return type treated as bool
	{
		return s >> (mer_length*2);
	}

	void get_neighbors(vector<k_mer> & collector)
	{
		for(int i = 0; i < value().ambi_nbs.size(); i++)
		{
			collector.push_back(get_neighbor(id, value().ambi_nbs[i].nb_info));
		}
		for(int i = 0; i< value().contig_nbs.size(); i++)
		{
			k_mer nid = value().contig_nbs[i].nid;
			if(nid != NULL_MER) collector.push_back(nid);
		}
	}

	void set_preds(MessageContainer & msgs) //assume ambiguous vertices have broadcasted msgs
	{
		//ambiguous vertices vote to halt; //unambiguous vertices set "preds" properly
		if(value().type == V_1)
		{
			k_mer self = (id | END_MER);
			value().preds.push_back(self); //one pred is itself
			//another is the only neighbor
			if(msgs.size() == 1)
			{
				value().preds.push_back(self);
			}
			else
			{
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(nbs.size() == 0)
					value().preds.push_back(self);
				else
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
				if(msgs[0] == nbs[0])
				{
					if(nbs.size() == 1)
						value().preds.push_back(self);
					else
						value().preds.push_back(nbs[1]);
				}
				else value().preds.push_back(nbs[0]);
			}
			else //msgs.size() == 0
			{
				k_mer self = (id | END_MER);
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(nbs.size() == 0)
				{
					value().preds.push_back(self);
					value().preds.push_back(self);
				}
				else if(nbs.size() == 1)
				{
					value().preds.push_back(self);
					value().preds.push_back(nbs[0]);
				}
				else
				{
					value().preds.push_back(nbs[0]);
					value().preds.push_back(nbs[1]);
				}
			}
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
		if(! Amb_is_SV)
		{
			if(step_num() == 1)
			{
				//ambiguous vertices send messages to neighbors
				if(value().type == Vm_n)
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
				if(value().type != Vm_n) LR_req(); //otherwise, there is no preds
			}
			//type 3 vertices will never receive any msg
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

class AmbLRAgg:public Aggregator<AmbLRVertex, bool, bool>
{
private:
	bool AND;
	long long int msgs_size;

public:
	AmbLRAgg(long long int msgs)
	{
		msgs_size = msgs;
		AND = true;
	}

	virtual void init()
	{
		AND=true;
	}

	virtual void stepPartial(AmbLRVertex* v)
	{
		if(Amb_is_SV)
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
		if(!Amb_is_SV)
		{
			if(step_num() % 2 == 1)
			{
				if(msgs_size != get_step_msg_num())
					msgs_size = get_step_msg_num();
				else
				{
					Amb_is_SV = true;
					global_step_num = 0;
				}
			}
		}
		return &AND;
	}

	virtual bool* finishFinal()
	{
		if(!Amb_is_SV)
		{
			if(step_num() % 2 == 1)
			{
				if(msgs_size != get_step_msg_num())
					msgs_size = get_step_msg_num();
				else
				{
					Amb_is_SV = true;
					global_step_num = 0;
				}
			}
		}
		return &AND;
	}
};

class AmbLRWorker:public Worker<AmbLRVertex, AmbLRAgg>
{
	char buf[100];

public:

	AmbLRWorker(int k)
	{
		set_mer_length(k);
		END_MER = 1ull << (2*mer_length);
	}

	virtual AmbLRVertex* toVertex(char* line)
	{
		char * pch;
		AmbLRVertex* v=new AmbLRVertex;
		pch=strtok(line, "\t");
		v->id = strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		v->value().type = (u8)atoi(pch);
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
		return v;
	}

	virtual void toline(AmbLRVertex* v, BufferedWriter& writer)
	{
		AmbLRValue & value = v->value();
		k_mer pred;
		if(value.type != Vm_n)
		{
			k_mer pred1 = (value.preds[0] & kick);
			k_mer pred2 = (value.preds[1] & kick);
			if(pred1 > pred2)
			{
				pred1 = pred2;
			}
			pred = pred1;
		}
		else
		{
			pred = NULL_MER;
		}
		int size = value.ambi_nbs.size();
		sprintf(buf, "%llu\t%llu %u %d", v->id, pred, value.type, size);
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
		writer.write(" #\n");
	}
};

void AmbListRank(string in_path, string out_path, int length)
{
	WorkerParams param;
	param.input_path=in_path;
	param.output_path = out_path;
	param.force_write=true;
	param.native_dispatcher=false;
//	init_workers();
	AmbLRWorker worker(length);
	AmbLRAgg agg = AmbLRAgg(-1);
	worker.setAggregator(&agg);
	worker.run(param);
//	worker_finalize();
}

#endif /* AMBILISTRANK_H_ */
