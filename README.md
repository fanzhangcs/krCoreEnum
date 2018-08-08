# krCoreEnum
The code is for the maximal (k,r)-core enumeration algorithm, published in the paper "When Engagement Meets Similarity: Efficient (k,r)-Core Computation on Social Networks", Fan Zhang, Ying Zhang, Lu Qin, Wenjie Zhang, Xuemin Lin, PVLDB 2017

# files
krCoreEnum.cpp - source code
data_edge.txt - toy friendship data with 1166 vertices and 1314 edges - data structure: vid \t nid \n...
data_attri.txt - vertex attribute data - geo-locations - data structure: vid \t latidude \t longitude \n...
//the data files are a part of the Gowalla dataset from SNAP: https://snap.stanford.edu/data/

# compile and run
complie with g++ and -O3
run and input the values of 'k' and 'r', such as 5 10 for k=5 and r=10km
the program ouputs in result.txt
