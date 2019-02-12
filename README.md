# findAntibodies
A custom Python program to accompany Meyer et al. "A Simplified Workflow for Monoclonal Antibody Sequencing." In submission. 2019.

## Description of the program: 
The Python program findAntibodies.py can be used as a guide to look for mouse antibody variable regions in the sequencing data. It is possible that some viable variable region sequences might not be properly identified by the program, especially if using the program to identify variable region sequences from organisms other than mouse. The program findAntibodies.py takes in a text file containing DNA sequencing data in FASTA format and identifies any kappa, lambda, and heavy variable regions in all frames of all the input DNA sequences by looking for common amino acid motifs found in these chains. The program returns which DNA sequences and which frame contain these motifs as well as which motifs were found and finally prints the full-length amino acid sequence of the found variable region. 
  
## Instructions for use  
  To use the program, put the Python file findAntibodies.py in the same folder as the text file with the DNA sequencing data in FASTA format. Rename the sequencing file DNA_seqs.txt. Then open a terminal and navigate to the folder containing the Python file and the text file using the cd (change directory) command. Once you have reached the proper folder, type:
  
#### findAntibodies.py <DNA_seqs.txt

and hit enter to view the results of the program in the command prompt window. Alternatively, you can type: 

#### findAntibodies.py <DNA_seqs.txt >Found_mAb_seqs.txt 

to view the results of the program in a new text file called Found_mAb_seqs.txt which will be created in the same working folder.
