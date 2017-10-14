#ifndef AMBISV_H_
#define AMBISV_H_

#include "utils/type.h"
#include "basic/DNAPregel-dev.h"
#include "GlobalDna.h"

using namespace std;

struct AmbiSVValue
{
	u8 type;
	k_mer prev_D;
	k_mer D;
	vector<k_mer> neighbors;
	vector<AmbiNB> ambi_nbs;
	vector<ContigNB> contig_nbs;
};

ibinstream & operator<<(ibinstream & m, const AmbiSVValue & v)
{
	m << v.type;
	m << v.prev_D;
	m << v.D;
	m << v.neighbors;
	m << v.ambi_nbs;
	m << v.contig_nbs;
	return m;
}

obinstream & operator>>(obinstream & m,  AmbiSVValue & v)
{
	m >> v.type;
	m >> v.prev_D;
	m >> v.D;
	m >> v.neighbors;
	m >> v.ambi_nbs;
	m >> v.contig_nbs;
	return m;
}

class AmbiSVVertex: public Vertex<k_mer, AmbiSVValue, k_mer>
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

	void set_neighbors(MessageContainer & msgs) //assume ambiguous vertices have broadcasted msgs
	{
		//ambiguous vertices vote to halt; //unambiguous vertices set "preds" properly
		if(value().type == V_1)
		{
			k_mer self = (id | END_MER);
			value().neighbors.push_back(self); //one pred is itself
			//another is the only neighbor
			if(msgs.size() == 1)
			{
				value().neighbors.push_back(self);
			}
			else
			{
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(nbs.size() == 0)
					value().neighbors.push_back(self);
				else
					value().neighbors.push_back(nbs[0]);
			}
		}
		else if(value().type == V1_1)
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
				if(msgs[0] == nbs[0])
				{
					if(nbs.size() == 1)
						value().neighbors.push_back(self);
					else
						value().neighbors.push_back(nbs[1]);
				}
				else value().neighbors.push_back(nbs[0]);
			}
			else //msgs.size() == 0
			{
				k_mer self = (id | END_MER);
				vector<k_mer> nbs;
				get_neighbors(nbs);
				if(nbs.size() == 0)
				{
					value().neighbors.push_back(self);
					value().neighbors.push_back(self);
				}
				else if(nbs.size() == 1)
				{
					value().neighbors.push_back(self);
					value().neighbors.push_back(nbs[0]);
				}
				else
				{
					value().neighbors.push_back(nbs[0]);
					value().neighbors.push_back(nbs[1]);
				}
			}
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
	}//in fact, a combiner with MIN operator can be used here

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
		value().D  =msgs[0];  //Once update the D[v], we should also update the Pre_D to keep the pre_D[v]
	}

	virtual void compute(MessageContainer & messages)
	{
		int cycle = 7;
		if(step_num() == 1)
		{
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
			set_neighbors(messages);
			if(value().type != Vm_n)
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
			if(value().type != Vm_n)
			{
				rtHook_2R(messages);
				rtHook_2S();
			}
		}
		else if(step_num()% cycle == 4)
		{
			if(value().type != Vm_n)
			{
				rtHook_3GDS(messages);
			}
		}
		else if(step_num()% cycle == 5)
		{
			if(value().type != Vm_n)
			{
				rtHook_4GD(messages);
			}
		}
		else if(step_num()% cycle == 6)
		{
			if(value().type != Vm_n)
			{
				rtHook_1S();
			}
		}
		else if(step_num()% cycle == 0)
		{
			if(value().type != Vm_n)
			{
				rtHook_2R(messages);
			}
		}
		else if(step_num()% cycle == 1)
		{
			if(value().type != Vm_n)
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

class AmbiSVAgg:public Aggregator<AmbiSVVertex, bool, bool>
{
private:
	bool AND;

public:
	AmbiSVAgg()
	{
		AND = true;
	}

	virtual void init()
	{
		AND=true;
	}

	virtual void stepPartial(AmbiSVVertex* v)
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

class AmbiSVWorker:public Worker<AmbiSVVertex, AmbiSVAgg>
{
	char buf[100];

public:

	AmbiSVWorker(int k)
	{
		set_mer_length(k);
		END_MER = 1ull << (2*mer_length);
	}

	virtual AmbiSVVertex* toVertex(char* line)
	{
		char * pch;
		AmbiSVVertex* v=new AmbiSVVertex;
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
		v->value().prev_D = v->id;
		v->value().D = v->id;
		return v;
	}

	virtual void toline(AmbiSVVertex* v, BufferedWriter& writer)
	{
		AmbiSVValue & value = v->value();
		k_mer pred;
		if(value.type != Vm_n)
		{
			pred = value.D;
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

void AmbiSV(string in_path, string out_path, int length)
{
	WorkerParams param;
	param.input_path=in_path;
	param.output_path = out_path;
	param.force_write=true;
	param.native_dispatcher=false;
//	init_workers();
	AmbiSVWorker worker(length);
	AmbiSVAgg agg = AmbiSVAgg();
	worker.setAggregator(&agg);
	worker.run(param);
//	worker_finalize();
}



#endif /* AMBISV_H_ */
