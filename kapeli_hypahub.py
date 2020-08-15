#!/usr/bin/python

import sys
import argparse

'''
HypaHub coding session problem

'''


def extract_kmers(sequence, k):
    '''
    Finds the set of kmers in a sequence
    '''
    subseq = []
    
    for ix,val in enumerate(sequence):
        
        if ix < (len(sequence)-k+1):
            subseq.append(sequence[ix:ix+k])
    
    return set(subseq)
           

def map_index_kmers(seq_list, k):
    '''
    Creates a dictionary of index to sequence (key) and the set of k-mers (value)
    '''
    k_dic = {}

    for idx, seq in enumerate(seq_list):

        k_dic[idx] = extract_kmers(seq, k) # get all kmers in sequence

    return k_dic


def reverse_complement(sequence):
    '''
    Returns the reverse complement of a sequence    
    '''
    comp = ""
    
    # substitute purine and pyrimidine pairs
    for s in sequence:
        if s == 'A':
            comp += 'T'
        elif s == 'T':
            comp += 'A'
        elif s == 'G':
            comp += 'C'
        elif s == 'C':
            comp += 'G'
        else:
            comp += s    # for 'N' or other non-ATCG characters
            
    # reverse the sequence
    comp = comp[::-1]
    
    return comp


def get_reverse_comp(seq_list):
    '''
    Returns a list of the reverse complement of sequences
    '''
    
    reverse_seq_list = []
    
    for x in seq_list:
        reverse_seq_list.append(reverse_complement(x))
       
    return reverse_seq_list


def print_dicts(seq, rev_comp):
    '''
    Print dictionaries of indexs mapped to k-mer set
    '''
    print("Sequence k-mers:")
    for k,v in seq.items():
        print(k, v)

        
    print("\nReverse complement k-mers:")
    for k,v in rev_comp.items():
        print(k, v)


def main(seq_list, k):
    '''
    Generates a list of integers represeting the index of sequences that contain at least one 
    k-mer present in another sequence (including the reverse complement).
    '''

    flagged = []
    
    # dictionary mapping kmers to index of sequences
    fwd_kmer_dic = map_index_kmers(seq_list, k)
    
    # dictionary mapping kmers to index of reverse complement of sequences
    reverse_seq_list = get_reverse_comp(seq_list)   
    rev_kmer_dic = map_index_kmers(reverse_seq_list, k)
    
    for idx, seq in enumerate(seq_list):

        for i in range(len(seq_list)):
            if idx != i:                                               # match with non-self sequences
                if fwd_kmer_dic[idx].intersection(fwd_kmer_dic[i]):    # if sequence has matching k-mer in another sequence
                    flagged.append(idx)                                # save index position of seqence and move to next sequence
                    break
                elif fwd_kmer_dic[idx].intersection(rev_kmer_dic[i]):  # if sequence has matching k-mer in reverse complement of another sequence
                    flagged.append(idx)                      
                    break

    print_dicts(fwd_kmer_dic, rev_kmer_dic)
    
    print("\nFlagged indexes:", flagged)



if __name__ == "__main__":
    
    # Parse data from file
    file = (sys.argv[1])
    
    with open(file) as f:
        lines = f.readlines() # list containing lines of file
    
    k = int(lines[0].strip('\n'))
    
    seq_list = []
    
    for i in range(1, len(lines)):
        s = lines[i].strip('\n')
        seq_list.append(s)
    
        
    # Check arguments and run main
    if type(k) is not int:
        print("Error: k is not of type int")
    else:
        main(seq_list, k)