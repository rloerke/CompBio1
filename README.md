# Computational Biology Assignment 1
## Written by Ray Loerke

Implementation of the Needleman-Wunsch global sequence alignment algorithm.
This program accepts a FASTA (.fasta) file containing the two sequences to be aligned.
It also uses a Matrix file (.mat) containing the scoring matrix to be used for the alignment.
If not Matrix file is provided default values will be used. 
The completed alignment will be written to an Alignment file (.g.out), which also displays statistics about the alignment.

See the example files for the formatting of these file types. 
Generally, a leading '<' means a sequence follows and a leading '#' means a comment follows.
