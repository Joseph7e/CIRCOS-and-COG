#!/usr/bin/python3
# Purpose: privided a list of parsed blast results. Assign COG terms using accesions to gi, gi to cog, and cog to name.
# Author Joseph Sevigny

import re

accession_gi = '/home/unhCW/xu/joe/cog_annotations/parsed_blast.txt'
gi_to_cog_file = '/home/unhCW/xu/joe/cog_annotations/cog2003-2014.csv'
cog_to_name = '/home/unhCW/xu/joe/cog_annotations/cog_names.tsv'
accession_to_loci_file = '/home/unhCW/xu/joe/cog_annotations/headers.txt'


cog_name_dict = {}
gi_cog_dict = {}
accession_gi_dict = {}
accession_locus_dict = {}


for line in open(accession_to_loci_file,'r').readlines():
    #>lcl|NZ_JWSS01000001.1_cds_WP_050926798.1_1 [gene=gltD] [locus_tag=RK51_RS00005] [protein=glutamate synthase] [fram
e=3]
    #print (line)
    line = line.rstrip()
    try:
        a = re.findall(r"(WP_[0-9]*).[1-9]_", line)[0]

    except:
        okay = "this wasn't a protein"
        #continue # this wasn't a protein

    l = re.findall(r"locus_tag=(.*)\]??", line)[0].split(']')[0]
    #print (a, l)
    accession_locus_dict[a] = l

for line in open(cog_to_name, 'r', encoding="latin-1").readlines():
    #COG0001 H       Glutamate-1-semialdehyde aminotransferase
    line = line.rstrip()
    if line[0] != '#':
        cog, letter, product = line.rstrip().split('\t')
        #print (cog,letter,product)
        cog_name_dict[cog] = [letter,product]

for line in open(accession_gi,'r').readlines():
    #WP_005495971.1 gi|15641976|ref|NP_231608.1| 296 264 263 245 57.55
    elements = line.rstrip().split(' ')
    accession = elements[0]; gi = elements[1].split('|')[1]
    if '.' in accession:
        accession = accession.split('.')[0]
    accession_gi_dict[accession] = gi

for line in open(gi_to_cog_file,'r').readlines():
    #158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
    elements = line.rstrip().split(',')
    gi = elements[0]; cog = elements[6]
    #print (gi, cog)
    gi_cog_dict[gi] = cog




count = 0

for a in accession_gi_dict.keys():
    try:
        accession = a
        cog = gi_cog_dict[accession_gi_dict[a]]
        letter = cog_name_dict[gi_cog_dict[accession_gi_dict[a]]][0]
        product = cog_name_dict[gi_cog_dict[accession_gi_dict[a]]][1].replace(' ','_')
        locus = accession_locus_dict[a]
        print (locus,cog,letter,product)
    except:
        cog = 'COG5411' # this one is not in the names dict for some reason
