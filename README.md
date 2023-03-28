# Stochastic_Resource_Control_to_AutoScale_MEC
This project introduces a stochastic resource control optimization model to allocate processing resources (CPUs/GPUs) at the edge for a case study of RAN slicing workload, i.e., LDPC decoding algorithm. LDPC decoding algorithm is modeled by asymtotic analysis to determine number of arithmetic operations. Peak performance of edge architecture is modeled by using Roofline parameters. CPLEX OPLIDE (.mod) and data (.dat) files are included.

**Pre-requisites**

--> g++ compiler > v7

--> jre > 1.5 for java compiler

--> IBM ILOG CPLEX (Educational License) ... You can obtain educational license from 
                                           https://www.ibm.com/academic/topic/data-science?ach_id=6fe17098-43df-4a9d-8412-3377286841a3


CPLEX was installed on x86 64-bit Ubuntu 20.04.4 LTS. Version for Windows and Mac machines are also available on IBM website.


>> Run stochastic model (joint_stoch.mod) with different values of datasets that are defined in the associated data file joint_stoch.dat). 

>> The stochastic model provides the results of S-RC that you can find in the Frontiers paper. 

>> Following that, run the the deterministic model (joint_determ.mod) using different values of R for each input data set as described in the paper. 

>> Values of data parameters are available in (joint_determ.dat). 
>> The output of running deterministric models are used to compute the probablistic processing time violation and overhead cost metrics as presented in the paper.   
