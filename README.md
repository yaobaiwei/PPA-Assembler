# PPA-Assembler  

De novo genome assembly is the process of stitching short DNA sequences to generate longer DNA sequences, without using any reference sequence for alignment.

PPA-assembler, a distributed toolkit for de novo genome assembly based on Pregel, a popular framework for large-scale graph processing. PPA-assembler adopts the de Bruijn graph based approach for sequencing and formulates a set of key operations in genome assembly. We implement these operations as **Practical Pregel Algorithms** (PPAs), which provide strong performance guarantees due to the bounds on computation and memory. The operations can also be flexibly assembled to implement various sequencing strategies according to users’ combination.  


## Highlights

* The first genome assembler based on Pregel, whose vertex-centric model is naturally fit for de novo genome assembly.
* We formulate a set of key operations that can be flexibly assembled to implement various sequencing strategies, where each genome assembly operation is a **PPA**, which provides strong performance guarantee.
* PPA-assembler demonstrates obvious advantages on efficiency, scalability, and sequence quality, comparing with existing distributed assemblers (e.g., **ABySS**, **Ray**, **SWAP-Assembler**).


## Getting Started

* **Install**  
  PPA-assembler is built on the top of our previous project [Pregel+](http://www.cse.cuhk.edu.hk/pregelplus/index.html). To install PPA-assembler's dependencies (e.g., MPI, HDFS), using instructions in the following [guide](http://www.cse.cuhk.edu.hk/pregelplus/documentation.html).

* **Build**   
	```bash
	$cd ${PPA_ROOT}/
	$./auto-build.sh
	```

* **Run**  
	```bash
	$cd ${PPA_ROOT}/release
	$mpiexec -n 1 ./put INPUT_FASTQ_FILE_PATH OUTPUT_HDFS_PATH
	$mpiexec -f /path/to/machine.conf -n M ./run
	```

* [Tutorials](docs/TUTORIALS.md)


## Academic and Reference Papers

[**VLDB 2014**] [Pregel Algorithms for Graph Connectivity Problems with Performance Guarantees](docs/ppa-vldb2014.pdf). Da Yan, James Cheng, Kai Xing, Yi Lu, Wilfred Ng, Yingyi Bu. PVLDB, Volume 7(14), Pages 1821-1832.

[**ICDE 2018**] Scalable De Novo Genome Assembly Using Pregel. Da Yan, Hongzhi Chen, James Cheng, Zhenkun Cai, Bin Shao. In Proceedings of the 34nd IEEE International Conference on Data Engineering (2018). [Full Paper Version](docs/ppa-assembler.pdf)
