# KLP
## Compile & Run

The source code of our paper is in the folder "src".

The source code of each algorithm for MalKLP enumeration (i.e., KLPE+, KLPE-BK+, KLPE-BR, KLPE-GR) has four parameters, corresponding to the dataset name (e.g., **D1**) and the three positive integers $k$, $l$ and $\theta$  (e.g., **(2,4,8)**).

* Compile

  > g++ main.cpp graph.cpp -O3 -o main


* Run
  
  > ./main D1 2 4 8


The source code of each algorithm for MaxKLP identification (i.e., MKLP+, MKLP-BS, MKLP-GR) has three parameters, corresponding to the dataset name (e.g., **D1**) and the two positive integers $k$ and $l$  (e.g., **(2,4)**).

* Compile

  > g++ main.cpp graph.cpp -O3 -o main


* Run
  
  > ./main D1 2 4


## Datasets

We provide five sample datasets in the folder "sample_dataset" for testing, D1 - D5.
