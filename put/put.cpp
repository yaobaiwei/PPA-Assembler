#include "utils/hdfs_core.h"

int main(int argc, char** argv)
{
	const char* input = argv[1];  	//input path on local, a fastq file
	const char* output = argv[2]; 	//output path on HDFS, as input of PPA-Assembler
	putFASTQ(input, output);
	return 0;
}
