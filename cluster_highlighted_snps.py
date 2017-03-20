#!/usr/bin/python3
# Purpose: Group highlighted snps into cluster groups
# Author: Joseph Sevigny

import sys

input = sys.argv[1] # output from parsing snp data
input_handle = open(input,'r')

gap_number = 10000 # number of nucleotides between snps to initiate new cluster

old_clustering_number = 0 # keep track of current clustering group
current_cluster = 1

for line in input_handle:
    line = line.rstrip()
    elements = line.split(',')
    total, red, yellow, id1, id2, product_name, chr, start, stop = line.split(',')
    start = int(start); stop = int(stop)
    current_clustering_number = start
    if current_clustering_number <= old_clustering_number+gap_number:
        print (line+','+str(current_cluster))
    else:
        print()
        current_cluster+=1
        print (line+','+str(current_cluster))
    old_clustering_number = stop
