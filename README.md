# deBruijn-Graph-Assembler
Python program to assemble a fragmented sequence using the deBruijn graph approach.

## Quick Start : 

#### Import the package : 
#### from deBruign_graph_assembler.py import *


db=deBruijnGraph(list) #e.g. list=['atgc', 'atgc', 'tcga', 'tgca', 'atcg']
###### OR
db=deBruijnGraph()
db.load_seq(sequence, k) #k is an integer for length of k-mers
db.assemble()


## This program contains three classes : DNA, ShortestCommonSuperstring and deBruijnGraph


### 1. DNA Class :
This class provides the functionality to take a sequence and break it down into k-mers, it can also provide unique k-mers.

#### 1.1 Initialising class and loading a sequence : 
dna=DNA(sequence)

#### 1.2 Finding all k-mers in the sequence (user has to provide length of k-mers = k) : 
all_kmers=dna.all_kmer(k) #k is an integer for length of k-mers

#### 1.3 Finding unique k-mers in the sequence : 
unique_kmers=dna.unique_kmers() #all_kmers() must be run in prior


### 2. ShortestCommonSuperstring Class (This class inherits from DNA class):
This class provides the functionality to make a directed graph by calculating edge weights from the overlaps of k-mers with each others. It can recursively and greedily merge k-mers with maximum overlap and reduce the list until either the SCS is found ot there are no further overlaps.

#### 2.1 Initialising class and loading a kmers (user has to provide list of k-mers) : 
scs=ShortestCommonSuperstring(list) #e.g. list=['ATGC', 'TGCC', 'GCCA']
scs is a list object that contains the shortest common superstring. It could have multiple strings if program is not able to resolve.

#### 2.2 Initialising class if user does not have a list of kmers (user has to provide a sequence) : 
obj1=ShortestCommonSuperstring()
You will get a message : Warning! No kmers provided. You can load sequences using load_seq() function.
obj1.load_seq(sequence, k) #k is an integer for length of k-mers

Finding SCS : 
scs=obj1.scs()


### 3. deBruijnGraph Class (This class inherits from ShortestCommonSuperstring class):
This class provides the functionality to make a dBruijn graph by calculating edge weights from the overlaps of k-1 mers with each others. It can recursively traverse across the Eulerian walk with maximum overlap and reduce the list until either the assembly is found or there are no further overlaps.

#### 3.1 Initialising class and loading a kmers (user has to provide list of k-mers) : 
db=deBruijnGraph(list) #e.g. list=['ATGC', 'TGCC', 'GCCA']
db.kmers is a list object that contains the final assembly. It could have multiple strings if program is not able to resolve.

#### 3.2 Initialising class if user does not have a list of kmers (user has to provide a sequence) : 
obj1=deBruijnGraph()
You will get a message : Warning! No kmers provided. You can load sequences using load_seq() function.
obj1.load_seq(sequence, k) #k is an integer for length of k-mers

Finding assembly : 
assembly=db.assemble()
