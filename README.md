MULKSG has the same requirements as SPAdes, but it also requires an MPI implementation (openmpi or mpich or similar) as well as mpi4py. To install mpi4py, 'pip install mpi4py'. You can also use conda to install the same way. Currently, uses the same options as SPAdes, you can use all the same inputs or yaml files. To run you just use mulksg.py instead of spades.py:  '/usr/bin/time -v mpirun -n 7 mulksg.py -1 ~/TestDataSets/Staphylococcus_aureus/raw/frag_1.fastq -2 ~/TestDataSets/Staphylococcus_aureus/raw/frag_2.fastq -o mulksg_Staphylococcus_aureus --careful -k 21,33,45,57,69,81,93'

Currently MULKSG will create all graphs, except the highest value k graph, at the same time. They all do error correction on the graphs and traverse to get contigs. Then the highest value k graph is built with the contigs generated. Currently the difference is that each graph is treated as if it were the last_one in the spades pipeline (This does extra error correction). The results are promising, as we can now do parallel computing. On one node we get about a 12% increase in core usage, with a small speedup(If there are sufficient resources on the node). But, we can now also assemble on multiple nodes. At this point I have not parallelized the pre-read correction or the post-contig correction, just the graph construction, traversal, and graph error correction. To make even a better tool, it would be best to complete parallelization for all stages.

As can be seen from the paper to AlCoB 2019 (https://drive.google.com/file/d/1VhTjF9L_5BRb_mtb6goHMAZ3jyUVGrKm/view), if there are limited resources on your desktop computer the speed is decreased. 

If parallel doesn't use the right number of threads, export OMP_NUM_THREADS XX, then add -x OMP_NUM_THREADS to the command line.


To build, cd to the spades_src folder, then type make. This will make both spades version 3.13.0 and mulksg.py. They share the same base code at this point.
