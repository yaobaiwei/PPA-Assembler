//#define SV_USED

#ifdef SV_USED
#include "dna/SV.h"
#include "dna/AmbiSV.h"
#else
#include "dna/ListRank.h"
#include "dna/AmbiListRank.h"
#endif

#include "dna/DeBruijn.h"
#include "dna/ContigMerge.h"
#include "dna/BubbleFilter.h"
#include "dna/AmbiConnect.h"
#include "dna/TipRemoval.h"
#include "dna/AmbiMerge.h"

int k_mer_t;
int freq_t;
int bubble_t;
int tip_t;
int output_contig_t;

string HDFS_INPUT_PATH;
string DeBruijn_PATH;
string KmerLink_PATH;
string AmbVtx_PATH;
string NoBubble_PATH;
string Contig_PATH;
string AmbConnect_PATH;
string NoTip_PATH;
string AmbLink_PATH;
string HDFS_OUTPUT_PATH;

void load_system_parameters()
{
	dictionary *ini;
	double val, val_not_found = -1;
	char *str, *str_not_found = "null";

	char* PPA_Assembler_HOME = getenv("PPA_Assembler_HOME");
	if(PPA_Assembler_HOME == NULL)
	{
		PPA_Assembler_HOME = "../";
	}
	string conf_path(PPA_Assembler_HOME);
	conf_path.append("./PPA_Assembler_conf.ini");
	ini = iniparser_load(conf_path.c_str());
	if(ini == NULL)
	{
		fprintf(stderr, "can not open %s!\nExits.\n", "PPA_Assembler_conf.ini");
		exit(-1);
	}

	// [PPA_Assembler]
	val = iniparser_getint(ini, "PPA_Assembler:k_mer_t", val_not_found);
	if(val!=val_not_found) k_mer_t=val;
	val = iniparser_getint(ini, "PPA_Assembler:freq_t", val_not_found);
	if(val!=val_not_found) freq_t=val;
	val = iniparser_getint(ini, "PPA_Assembler:bubble_t", val_not_found);
	if(val!=val_not_found) bubble_t=val;
	val = iniparser_getint(ini, "PPA_Assembler:tip_t", val_not_found);
	if(val!=val_not_found) tip_t=val;
	val = iniparser_getint(ini, "PPA_Assembler:output_contig_t", val_not_found);
	if(val!=val_not_found) output_contig_t=val;

	str = iniparser_getstring(ini,"PPA_Assembler:HDFS_INPUT_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) HDFS_INPUT_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:DeBruijn_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) DeBruijn_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:KmerLink_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) KmerLink_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:AmbVtx_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) AmbVtx_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:NoBubble_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) NoBubble_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:Contig_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) Contig_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:AmbConnect_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) AmbConnect_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:NoTip_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) NoTip_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:AmbLink_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) AmbLink_PATH = str;
	str = iniparser_getstring(ini,"PPA_Assembler:HDFS_OUTPUT_PATH", str_not_found);
	if(strcmp(str, str_not_found)!=0) HDFS_OUTPUT_PATH = str;

	iniparser_freedict(ini);
}

int main(int argc, char** argv)
{
	init_workers();
	load_system_parameters();

	//sample
	DeBruijn_Build(HDFS_INPUT_PATH, DeBruijn_PATH, k_mer_t, freq_t);  //freq threshold
	worker_barrier();

#ifdef SV_USED
	SV(DeBruijn_PATH, KmerLink_PATH, AmbVtx_PATH, k_mer_t);
	worker_barrier();

	Contig_Merge(KmerLink_PATH, NoBubble_PATH, k_mer_t);
	worker_barrier();

	Bubble_Filter(NoBubble_PATH,  Contig_PATH, k_mer_t, bubble_t);  //editDistance
	worker_barrier();

	AmbiConnect(AmbVtx_PATH, Contig_PATH,AmbConnect_PATH, k_mer_t, tip_t);  //tip's sequence length
	worker_barrier();

	TipRemoval(AmbConnect_PATH, NoTip_PATH, k_mer_t, tip_t);  //tip's sequence length
	worker_barrier();

	AmbiSV(NoTip_PATH, AmbLink_PATH, k_mer_t);
	worker_barrier();

#else
	ListRank(DeBruijn_PATH, KmerLink_PATH, AmbVtx_PATH, k_mer_t);
	worker_barrier();

	Contig_Merge(KmerLink_PATH, NoBubble_PATH, k_mer_t);
	worker_barrier();

	Bubble_Filter(NoBubble_PATH,  Contig_PATH, k_mer_t, bubble_t);  //editDistance
	worker_barrier();

	AmbiConnect(AmbVtx_PATH, Contig_PATH,AmbConnect_PATH, k_mer_t, tip_t);  //tip's sequence length
	worker_barrier();

	TipRemoval(AmbConnect_PATH, NoTip_PATH, k_mer_t, tip_t);  //tip's sequence length
	worker_barrier();

	AmbListRank(NoTip_PATH, AmbLink_PATH, k_mer_t);
	worker_barrier();

#endif
	AmbiMerge(Contig_PATH, AmbLink_PATH, HDFS_OUTPUT_PATH, k_mer_t, output_contig_t);

	worker_finalize();
	return 0;
}
