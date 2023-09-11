This zip includes:

- Source code of our 2 streaming algorithms, FSA and IFSA.
- Facebook and Sensor dataset (in the "data" folder) for testing the algorithms.

Since estimating F in influence maximization is extremely time-consuming, our code uses OpenMP for parallelization (https://en.wikipedia.org/wiki/OpenMP).

To build our code, run:

```
	g++ -std=c++11 *.cpp -o ksub -DIL_STD -fopenmp -g
```

After building, to run our code, run:

```
	./kSMK -f <data filename> 
		-V <size of V> 
		-t <type of experiment, 0: influence maximization, 1: sensor placement> 
		-k <value of k> 
		-B <value of B> 
		-e <value of epsilon>
		-n <value of eta - denoise step for RStream> 
		-a <algorithm, 0: FSA, 1: IFSA Please use SSA source code for testing IM algorithm>
		-p <number of threads (OpenMP) to running algorithms>

```

All the implementations are on Linux machine with configurations of 2 x Intel(R) Xeon(R) CPU E5-2697 v4 @ 2.30GHz and 10 threads x 128 GB DIMM ECC DDR4 @ 2400MHz. 