#!/usr/bin/env python3
#written by Lena Meyer
#February 2019 final version

'''Find antibody kappa, lambda and heavy variable region sequences
in Sanger sequencing data.

Input: A FASTA file containing DNA sequences
resulting from Sanger sequencing of a blunt-ended cloning product.
Output: A report on whether an antibody variable region sequence
was found in each input DNA sequence. It will also print the motifs found.'''

import sys
import re
import itertools

class FastaReader:
    '''Parse the input FASTA file.'''

    def __init__(self, fastaFile):
        '''fastaFile is a file object input by the user.'''
        
        self.fastaFile = fastaFile
        
    def readFasta(self):
        '''Yields FASTA file entries as (header, sequence).'''

        #initialize empty seq string
        seq = ''
        
        for line in self.fastaFile:
            #remove \n character from the end of each line
            line = line.strip()

            #lines starting with > are headers
            if line.startswith('>'):
                #if seq isn't empty, head + seq are both occupied,
                #ie. there is a FASTA entry to return, so
                #yield (head, seq)
                if seq != '':
                    yield (head, seq)
                    #retinitalize empty seq string after yielding
                    seq = ''
                #remove > from header    
                head = line[1:]
                
            #lines with alphabet letters are sequences    
            elif line.isalpha():
                #+= means seq = seq + line
                seq += line
                
        #make sure to yield the last entry
        #there is no > at the end of the FASTA file
        #so the last entry wouldn't be yielded otherwise
        yield (head, seq)
        

class MotifCompiler:
    '''Make lists of amino acid motifs found in kappa, lambda, and heavy chains.
Compile the motifs (means saving a regex pattern into a regex object) once here
so that compilation doesn't have to happen again for every sequence.'''
    
    def __init__(self):
        #make lists of amino acid motifs
        #found in the beginning, middle, and end (ie. J-region start)
        #of kappa, lambda, and heavy chain variable regions
        #represent motifs using regular expression notation
        self.kappaMotifs = ['[I|T|V][A-Z][L|M][T|S]Q[S|T|P][P|H|T|S]',
                            'P[S|D|A|V][R|H]F[S|T|R]GS[G|D|N|R]',
                            'TFG[G|Q|A|S|T]?GT']
        self.lambdaMotifs = ['[V|L][T|H]Q[E|S|P][S|P|A][A|S|L][L|A|V][T|S][T|S|F|G]',
                             '[P|D|S][G|D][V|I|L][P|S][A|V|D|P|N]RFSGS[L|K|S]',
                             '[F|W][G|E][G|S|T][G|E|R][T|K][K|R|S][L|V|Y|T][T|L][V|D|W]']
        self.heavyMotifs = ['V[Q|K|N|M]L[V|K|Q|L|H][Q|E][S|P]G',
                            'W[V|I][R|K][Q|K][A-Z][P|N|H]',
                            '[Q|E|A|H]G[T|S][L|T|S|M][V|L][T|A]VS[S|A]']

    def compileMotifs(self):
        '''Make a list, self.compiledKappaMotifs, with format [[list], 'single element'].
The first element will be a list of compiled kappa chain motifs (ie. patterns).
The second element will be the chain type string (ex. 'kappa') to keep track.
Do the same for lambda and heavy chain motifs.
Put all together in a list of list of lists called self.allCompiledMotifs,
which will allow for searching for all motifs at once within one sequence
and printing the result of the search.'''
        
        self.compiledKappaMotifs = [[re.compile(motif)
                                     for motif in self.kappaMotifs],
                                    'kappa']
        self.compiledLambdaMotifs = [[re.compile(motif)
                                     for motif in self.lambdaMotifs],
                                    'lambda']
        self.compiledHeavyMotifs = [[re.compile(motif)
                                     for motif in self.heavyMotifs],
                                    'heavy']
        self.allCompiledMotifs = [self.compiledKappaMotifs,
                                  self.compiledLambdaMotifs,
                                  self.compiledHeavyMotifs]

        return self.allCompiledMotifs

class FindAntibodies:
    '''Define a new class called FindAntibodies which will identify
possible kappa, lambda and heavy chains in all frames of input DNA sequences.'''

    #make a dictionary of {RNA codons: amino acids} and the DNA equivalent
    #DNA codon table made from RNA codon table by replacing the Us with Ts
    rnaCodonTable = {
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U','T'): value
                      for key, value in rnaCodonTable.items()}
    
    
    def __init__(self, head, seq, motifs):

        #accept and store sequences (as capital letters) and headers
        #also store the motifs compiled in class MotifCompiler
        self.sequence = seq.upper()
        self.header = head
        self.allCompiledMotifs = motifs
                
        #store the results of the translateSeq and getRevComp methods
        #ie. the AA translations of the six frames per sequence
        #in one list with six elements, self.fwdrevTranslations
        #will use this list in the whereAntibodies method
        #to search for antibody variable chain motifs
        self.fwdrevTranslations = self.translateSeq(self.sequence) + self.getRevComp()
       
                
    def translateSeq(self, seq):
        '''Translate DNA sequences assuming a forward 5' to 3' sequence.'''
        
        #make a list of 3 lists, fwdAAs, which will each contain
        #the amino acids of a frame as list elements (many elements)
        fwdAAs = [[], [], []]
        
        #initialize list joinedAAseqs which will contain
        #the joined translated DNA sequences of 3 frames (3 total elements)
        joinedAAseqs = []
            
        #loop through the forward frames and translate each frame
        for frame in range(0, 3):
            #go in steps of three to look at each codon
            for i in range(frame, len(seq), 3):
                #translate unknown amino acids (containing N in codon) as X
                #and append amino acids to fwdAAs[frame] to keep track of
                #which amino acids belong to which frame
                if 'N' in seq[i: i+3]:
                    thisAA = 'X'
                    fwdAAs[frame].append(thisAA)
                #translate known amino acids using the DNA codon dictionary and
                #append amino acids to fwdAAs[frame]
                elif seq[i: i+3] in self.dnaCodonTable.keys():
                    thisAA = self.dnaCodonTable.get(seq[i: i+3])
                    fwdAAs[frame].append(thisAA)
            #join the list elements in fwdAAs to make complete AA sequences
            #append this joined sequence to joinedAAseqs
            joinedAAseqs.append(''.join(fwdAAs[frame]))
        
        return(joinedAAseqs)
                          

    def getRevComp(self):
        '''Calculate the reverse complement of the entered sequence
and give the reverse complemented sequence to the translateSeq method.'''

        #make a dictionary of base pairs
        #unknown base N will be paired with itself
        bpDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}

        #get the value of the key for a base, put it in a list,
        #and join the list elements into a string
        #results in the reverse complement of the orignal sequence
        revComp = ''.join([bpDict[base] for base in reversed(self.sequence)])
       
        #take the reverse complemented DNA sequence and
        #translate to AA sequences using the translateSeq method
        revTranslations = self.translateSeq(revComp)
        return(revTranslations)


    def whereAntibodies(self, head, seq):
        '''Look for amino acid motifs of kappa, lambda, and heavy chains in
all six frames of a sequence and report if there are any
full-length or partial antibody variable region sequences.'''               

        #Overview:
        #for each translation found in the six frames of a DNA sequence,
            #look for AA motifs of all chain types
                #print information on whether any or all motifs were found
                #ie. whether a partial or a full-length chain was found
                #print motifs, >self.header, and the translation on a new line
                                       
        #search for all AA motifs of kappa, lambda, and heavy chains
        #motifs are a list given in the first element [0]
        #of each list in self.allCompiledMotifs
        #search in given DNA sequence translations for these motifs
        for translation in self.fwdrevTranslations:
            #the frame corresponds to the index of the translation
            #in the list self.fwdrevTranslations
            frame = self.fwdrevTranslations.index(translation)
           
            for motifs, chaintype in self.allCompiledMotifs:
                motifMatches = []
                for motif in motifs:
                    #append all motifs that are found in the translation to motifMatches
                    #motifMatches then looks like eg. [['motif'], ['motif'], []]
                    #if the first 2 of 3 motifs for a chain were found
                    motifMatches.append(motif.findall(translation))

                #create a new list, flattenedMotifMatches
                #to flatten all lists in motifMatches
                #ie. the ['motif'] lists and empty lists (if no motif was found)
                #itertools.chain(*iterables) is the general format
                flattenedMotifMatches = itertools.chain(*motifMatches)
               
                #if all sequence motifs are found,
                #ie. each list in motifMatches is occupied with a motif string,
                #the sequence is a full-length antibody variable region sequence
                if all(motifMatches):                    
                    print("Complete {} chain found in frame "
                              .format(chaintype), end='')
                    #for forward frames
                    if frame in range(0, 3):
                        print("{} of {}! \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        frame+1))
                    #for reverse frames
                    else:
                        print("{} of {}! \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        -frame+2))
                    print("{}\n".format(translation))
                #if any sequence motifs are found,
                #ie. some but not all lists in motifMatches are occupied,
                #the sequence is a partial antibody variable region sequence
                elif any(motifMatches):                                               
                     print("Partial {} chain found in frame "
                              .format(chaintype), end='')
                    #for forward frames
                     if frame in range(0, 3):
                        print("{} of {}. \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(frame+1,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        frame+1))
                    #for reverse frames
                     else:
                        print("{} of {}. \nMotifs found: \n{}\nAmino acid sequence in FASTA format: \n>{}_frame_{}_translation"
                                .format(-frame+2,
                                        self.header,
                                        ', '.join(flattenedMotifMatches),
                                        self.header,
                                        -frame+2))
                     print("{}\n".format(translation))
                    
                                                                 
#instantiate the MotifCompiler object
motifs = MotifCompiler().compileMotifs()

#read the input FASTA file and instantiate the FindAntibodies object
#using the FASTA file info and the MotifCompiler object
myReader = FastaReader(sys.stdin)
for head, seq in myReader.readFasta():
    checkAllSequences = FindAntibodies(head, seq, motifs)
    checkAllSequences.whereAntibodies(head, seq)
