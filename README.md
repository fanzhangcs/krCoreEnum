# krCoreEnum
The code is for the maximal (k,r)-core enumeration algorithm, published in the paper "When Engagement Meets Similarity: Efficient (k,r)-Core Computation on Social Networks", Fan Zhang, Ying Zhang, Lu Qin, Wenjie Zhang, Xuemin Lin, PVLDB 2017

# files
krCoreEnum.cpp - source code 

data_edge.txt - toy friendship data - data structure: vid \t nid \n... - note that each edge is stored twice and ordered here

data_attri.txt - vertex attribute data - geo-locations - data structure: vid \t latidude \t longitude \n...

the data files are a part of the Gowalla dataset from SNAP: https://snap.stanford.edu/data/


# compile and run
complie with g++ and -O3

run and input the values of 'k' and 'r', such as 5 10 for k=5 and r=10km

the program ouputs in result.txt

# note
If you have any question, please contact me by fanzhang.cs@gmail.com.

If you used this code, please kindly cite the paper.
