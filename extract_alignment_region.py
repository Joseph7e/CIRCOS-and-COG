#!/usr/bin/python3
# Purpose: Given a set of coordinates, extract portions of an alignment.
# Author: Joseph Sevigny


import sys
from Bio import SeqIO

out_add = sys.argv[3]
list_starts = []
list_stops = []


for line in open(sys.argv[1],'r').readlines():
    number, start, stop = line.rstrip().split(' ')
    list_starts.append(start)
    list_stops.append(stop)

print (list_starts)
print (list_stops)


ref_id = '10329'
alignment_file = sys.argv[2]

reference = ''
non_references = {}



for seq_record in SeqIO.parse(alignment_file, "fasta"):
    if ref_id in str(seq_record.id):
        reference = str(seq_record.seq)
    else:
        current_header = str(seq_record.id)
        current_seq = str(seq_record.seq)
        non_references[current_header] = current_seq

alignment_location = -1
real_location = 0

alignment_starts = []
alignment_stops = []

for nucleotide in reference:
    alignment_location += 1
    if not nucleotide == '-':
        real_location += 1
    if str(real_location) in list_starts:
        alignment_starts.append(alignment_location)
    if str(real_location) in list_stops:
        alignment_stops.append(alignment_location)

print (alignment_starts)
print (alignment_stops)

for i in range(len(alignment_starts)):
    outfile = open('partial_alignment_'+out_add+'_'+str(i+1)+'.fasta','w')
    aln_start = alignment_starts[i] ; aln_stop = alignment_stops[i]
    real_start = list_starts[i]; real_stop = list_stops[i]
    outfile.writelines('>'+ref_id+'_aln_'+str(aln_start)+'_'+str(aln_stop)+'\n'+reference[aln_start:aln_stop]+'\n')
    for head, seq in non_references.items():
        outfile.writelines('>'+head+'_aln_'+str(aln_start)+'_'+str(aln_stop)+'\n'+seq[aln_start:aln_stop]+'\n')
