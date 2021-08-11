# Vectorizing BigSMILES into a tuples counter vector

## What are the tuples

The tuples are all the singlets, doublets and triplets that can be formed by a given list of non-hydrogen atomes.
A singlet is a non-hydrogen atom is characterized by the atom type, hybridization (either 'sp3', 'sp2' or 'sp1'), the number of non-hydrogen atoms it is bonded to, and its position in the repeat unit (backbone or side chain). 
A doublet (triplet) consists of two (three) singlets that have a covalent bond between them specified by one of the following connectors: ‘’, ‘=’ or ‘:’. The singlets, doublets and triplets are then concatenated into one large feature vector.

## Generation of tuple counter vectors

For each polymer in the dataset, the singlets, doublets, and triplets are generated using an infinitely long polymer representation, encoded using the BigSMILES representation. The polymer's repeat unit is made into a graph representation, including the connector type (‘$’ or ‘<’ or ‘>’). From this graph representation, all the singlets are found by looking at each node of the graph. Similarly for doublets and triplets, the nodes and edges of the graph are scanned through to determine the bond type between atoms. The tuples vector starts as a vector of zeros of the length of the feature vector. For each tuples that occurs in the polymer, a count is added to the corresponding tuple's index. The final tuples vector of singlets, doublets and triplets is the occurrence count of each tuple divided by the total number of non-hydrogen atoms in a repeat unit of the polymer. 
Once all the polymers have been vectorized in this way, a tuples vector matrix is generated with each column corresponding to a polymer of the dataset and each row corresponding to a tuple. The last step is the reduction of the matrix by removing any row that consists of only zeros (meaning this tuple did not occur in any of the polymers of the dataset). A corresponding string format list of the tuples is generated, which keeps the same ordering as the tuples matrix.

## Input and output of the code

The code to run on Python3 to generate the tuple vecotr is BigSMILES_vectorization.py.

### Code input

A CSV file with a single column containing the BigSMILES string of the polymer you wish to vectorize. 
Be sure to write each BigSMILES string in the correct format. For examples, polypropylene would be written as: '{[][$]CC(C)[$][]}'.

### Code output

The output of the code is 6 files.
The full tuple counter vector of the polymer dataset will be saved as a CSV file, with a file for the reduced singlets, doublets and triplets respectively.
For each of the tuple vector files, a CSV file with the correctly ordered string tuples will be generated. For example, the propylene C-C bond containing the side group carbon will be written as: [['C', 'sp3', 3, 1, 'sg'], '-', ['C', 'sp3', 1, 3, 'bb']], where the singlet ['C', 'sp3', 3, 1, 'sg'] specifies, in order, the atom type, the hybridization, the number of neighboring hydrogen atoms, the number of neighboring non-hydrogen atoms, and the position in the repeat unit.

### Running the code

In the python script BigSMILES_vectorization.py, use the function vectorize that takes two inputs: the csv file and the location of that file in your system. (see example at end of script)

### Code limitations

Currently, the only atoms that can be present in the polymer dataset are C, c, N, n, O, o, Si, S, F, Cl, Br (none capitalized being aromatic atoms). The code should only be used for linear polymers. Repeat units that have more than one connector would lead to inproper counting of the bonds that connect different repeat units. The code is currectly only effective for polymers with a single repeat unit. This vectorization scheam cannot account for specific stereoisomerism.




