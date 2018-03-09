#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 16:36:28 2018

@author: anton
"""
from Bio import SeqIO

class Kmer:
    counter = 0
    sequence = ''
    
    def __init__(self, kmer_name):
        self.sequence = kmer_name
        self.position_list = []
    
    def increase(self):
        self.counter += 1
    
    def addpos(self, chrom, start, end):
        self.position_list.append((chrom, start, end))
    
    def statistics(self):
        print('statistics for ',self.sequence,' Kmer')
        print ('chromosome', '\t', 'cor1', '\t', 'cor2', '\t')
        for locus in self.position_list:
            print (locus[0], '\t', locus[1], '\t', locus[2])


kmer_size = 23
kmer_dict = {}

for record in SeqIO.parse('/home/anton/seq_y_pestis.fasta', 'fasta'):
    seq = str(record.seq)
    seq_name = str(record.id)
    seq_leng = len(record.seq)
    for index in range (seq_leng-kmer_size+1):
        print(index)
        current_kmer = seq[index:(index+kmer_size)]
        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].increase()
            kmer_dict[current_kmer].addpos(seq_name, index,index+kmer_size-1)
        else:
            kmer_dict[current_kmer] = Kmer(current_kmer)
            kmer_dict[current_kmer].increase()
            kmer_dict[current_kmer].addpos(seq_name, index,index+kmer_size-1)

maxmer = 0
for kmer in kmer_dict:
    if maxmer < kmer_dict[kmer].counter:
        maxmer = kmer_dict[kmer].counter
        maximer = kmer_dict[kmer]
        print(maximer.sequence)    
maximer.statistics()


