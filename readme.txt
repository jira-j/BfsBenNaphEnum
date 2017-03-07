BfsBenNaphEnum Manual

BfsBenNaphEnum is a fast enumeration tool that takes chemical formula in the form of the number of each atom types and enumerates chemical compounds without cycle except for benzene rings and naphthalene rings. 
You don't need to input the number of benzenen rings and naphthalene rings because this tool automatically calculates it for you. 

How to run 
Run "Make" command
Run "./bfsenum " command follows by the options that you want
    "-c" + number of carbon atoms
    "-n" + number of nitrogen atoms
    "-o" + number of oxygen atoms
    "-h" + number of hydrogen atoms
    "-t" to write the enumerated compounds in SMILES format to the file named "output.smi", otherwise it returns only the number of enumerated compounds on the screen

Or run with the following options if the input is a left-heavy chemical tree
    "-f" + the name of a file storing a left-heavy chemical tree
    "-t" to write the enumerated compounds in SMILES format to the file named "output.smi", otherwise it returns only the number of enumerated compounds on the screen

Example 
Run "./bfsenum -c6 -h6" if the chemical formula is C6H6 and do not need to know the structures.

Please refer to http://sunflower.kuicr.kyoto-u.ac.jp/~jira/bfsenum/index.html for the format and an example of the input file storing a left-heavy chemical tree.

Output explanation
This work enumerates each possible combination of the number of benzene rings and naphthalene rings per round.
Output of each round consists of two main parts. 
1. the number of nodes in the enumerated compounds and 
2. the accumulated number of enumerated structure until the current round. 

Example
Output of "./bfsenum -c6 -h6" is shown below.

===================
amount of naphthalene rings: 0
amount of benzene rings: 0
amount of C: 6
amount of N: 0
amount of O: 0
amount of H: 6
accumulated #enumerated structure until 1 round is 15
===================
amount of naphthalene rings: 0
amount of benzene rings: 1
amount of C: 0
amount of N: 0
amount of O: 0
amount of H: 6
accumulated #enumerated structure until 2 round is 16
======Result=======
Total #enumerated structures:16

Chemical compounds are enumerated in 2 rounds.
In the first round, 15 chemical compounds containing 0 benzene ring and 0 naphthalene ring are enumerated.
In the second round, 1 chemical compounds containing 1 benzene ring and 0 naphthalene ring are enumerated, so accumulated number of enumerated compounds is 16.
The number of all enumerated compounds is concluded again at the last line of the output.

Caution 
If a "double free allocation" error is detected during the run time, please use jemalloc to allocate the memory to avoid such error.
