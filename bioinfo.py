#!/usr/bin/env python
# Author: Samantha Obermiller sobermil@uoregon.edu
# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning


DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    my_list: list = []
    init_list(my_list)
    for x in range (0,101):
        #print(x)
        lst.append(value)
        #print(value)
    return lst

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """calculates the average quality score of the whole phred string"""
    sumscore=0
    for value in phred_score:
        sumscore = sumscore + convert_phred(value)
        averagescore = sumscore / len(phred_score)
        return(averagescore)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters"
    
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def oneline_fasta(sequence):
    '''convert fasta with two lines of sequence into one line of sequence'''
    sequence=""
    header=""
    with open("./one_line.fa", "w") as output:
        with open("./test.fa", "r") as fh:
            for line in fh:
                line=line.strip('\n')
                if ">" in line:
                    if header != "":
                        output.write(header+"\n")
                        output.write(sequence+"\n")
                        sequence=""
                    header=line
                else:
                    sequence+=line
            output.write(header+"\n")
            output.write(sequence+"\n")
            return ('success')
def reverse_complement(seq) -> str:
    "create reverse compliment of a given sequence of DNA"
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    sequence=list(seq)
    reverse_comp=str()
    for letter in reversed(sequence):
        reverse_comp += complements[letter]
    return reverse_comp
def flag_revcomp(strand):
    "determine if bitwise flag in sam file indicates read is reverse complement"
    if ((strand & 16) == 16):
        rev_comp = "reverse"
    else:
        rev_comp = "forward"
    return rev_comp


if __name__ == "__main__":
    # write tests for functions above

    #Covert_Phred Check
    """Check that convert_phred returns the correct value for several different inputs"""
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("convert_phred working")

    #Qual Score Check
    assert convert_phred("I")+convert_phred("@") == 71
    assert convert_phred("C")+convert_phred("2") == 51
    print("quality score working")


    #VALIDATE BASE SEQ
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    print("Passed DNA and RNA tests")

    #GC CONTENT
    assert gc_content("GCGCGC") == 1, "messed up calc when all GC"
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")

    #reverse compliment check
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("TACG") == "CGTA"