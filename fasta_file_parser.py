#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 09:27:41 2019

@author: rahul
"""

def read_FASTA(fname):
	"""(str) -> (list of tuples)
		Reads a fasta file and returns a list of tuples containing
		the header and the sequence
	"""
	with open(fname) as f:
		data = f.readlines()
		sequences = []
		sequence,sequence_name = '',''
		
		#Iterate through values in data
		for line in data:
			#Check for header
			if line[0] == '>':
				#Append seqeunce_name,sequence as a tuple to the list if sequence is not empty
				if sequence != '' : sequences.append((sequence_name,sequence))
				sequence = ''						#Reinitialize sequence  
				sequence_name = line[1:].split('|')[1]	#Assign sequence_name
			#Store sequence corresponding to header
			else:
				sequence += line.strip()
		#Append last seqeunce_name,sequence as a tuple to the list
		sequences.append((sequence_name,sequence))
	#Return list of tuples
	return sequences
f1 = open('mergedfile_output.tsv','w')
if __name__ == '__main__':
    compositiontmh = {}
    lentmh = 0 
    for seq_name, seq in read_FASTA("/home/rahul/Downloads/Seq2Feature/mergedfile"):
        print(seq_name,seq)
        f1.writelines(seq_name+'\t'+seq+'\n')