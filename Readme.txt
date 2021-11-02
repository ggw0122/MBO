
Monarch Butterfly Optimization (MBO)

Gai-Ge Wang
May 10, 2015

email:  gaigewang@163.com
       gaigewang@gmail.com

The files in this zip archive are MATLAB m-files that can be used to study Monarch Butterfly Optimization algorithm.

MBO is the method that we invented and wrote about in the following paper:
Gai-Ge Wang, Suash Deb, and Zhihua Cui, Monarch Butterfly Optimization.
Neural Computing and Applications, in press.
DOI: 10.1007/s00521-015-1923-y

Note: 
1) I do not make any effort to rewrite some general codes, while I reuse some codes according to Prof. Dan Simon. We should stand the shoulder of giant, therefore I have more time to focus on our method-MBO. In the following, I will provide the detailed description about all the codes. 
2) The MATLB R2015a is used when implementing our method. 
3) The C++ code of Monarch Butterfly Optimization (MBO) can be found in the web: https://github.com/ggw0122/Monarch-Butterfly-Optimization
4) As discussed in the paper, the MBO cannot find the best solution for each run. Our research team will improve it and distribute the codes in our future research. 

The MATLAB files can be used to reproduce the results in the paper, or to do your own experiments. The paper and the software are available at https://github.com/ggw0122/Monarch-Butterfly-Optimization. The software is freely available for any purposes (it is on the Internet, after all) although I would of course appreciate an acknowledgement if you use it as part of a paper or presentation.

The MATLAB files and their descriptions are as follows:

Ackley.m: 
This is the benchmark functions discussed in the paper. You can use it as template to write your own function if you are interested in testing or optimizing some other functions. This code is modified according to Dan Simon. The original one is available at http://academic.csuohio.edu/simond/bbo.

MBO_Generation_V2.m, MBO_FEs_V2.m:
Monarch Butterfly Optimization algorithm. The fixed generations (iterations) and fixed Function Evaluations (FEs) are considered as termination condition for MBO_Generation_V2.m and MBO_FEs_V2.m, respectively. It can be used to optimize some function by typing, for example, the following at the MATLAB prompt:
>> MBO_Generation_V2(@Ackley);
This command would run MBO_Generation_V2 on the Ackley function (which is codified in Ackley). 
>> MBO MBO_FEs_V2(@Ackley);
This command would run MBO_FEs_V2 on the Ackley function (which is codified in Ackley). 

Init.m: 
This contains various initialization settings for the optimization methods. You can edit this file to change the population size, the generation count limit, the problem dimension, the maximum Function Evaluations (FEs), and the percentage of population of any of the optimization methods that you want to run. This code is modified according to Dan Simon. The original one is available at http://academic.csuohio.edu/simond/bbo.

ClearDups.m: 
This is used by each optimization method to get rid of duplicate population members and replace them with randomly generated individuals. This code is modified according to Dan Simon. The original one is available at http://academic.csuohio.edu/simond/bbo.

ComputeAveCost.m: 
This is used by each optimization method to compute the average cost of the population and to count the number of legal (feasible) individuals. This code is the same as Dan Simon. The original one are available at http://academic.csuohio.edu/simond/bbo.

PopSort.m: 
This is used by each optimization method to sort population members from most fit to least fit. This code is the same with Dan Simon. The original one is available at http://academic.csuohio.edu/simond/bbo. 

Conclude1.m, Conclude2.m: 
They are concludes the processing of each optimization method. It does common processing like outputting results. Conclude1.m and Conclude2.m are used in MBO_Generation_V2.m and MBO_FEs_V2.m, respectively. They are modified according to Dan Simon. The original one is available at http://academic.csuohio.edu/simond/bbo.

I hope that this software is as interesting and useful to you as is to me. Feel free to contact me with any comments or questions.


