#!/usr/bin/python3

#Author: Joseph Sevigny
#Purpose: The purpose of the program is too parse snp variant data and
#           Construct input for circos

import sys, os

input_file_handle = open(sys.argv[1],'r')

###Order in excel sheet
#Reference, ST36(pacific), ST631(1G), ST631(2), ST631(1PQ), ST36(atlantic), ST636

###original patterns (no longer used)
# patterns1=['1121212'] #ST36, and 631(clade 2, MAVP-Q) ####
# patterns2 = ['1121211', '1212122'] #ST36, 631(MAVP-Q), and 636 ### Light Orange


###OUTPUT STUFF
output_dir = 'circos_data/'
#if not os.path.exists(output_dir): #creates a directory for files
os.mkdir(output_dir)
gene_info_handle = open(output_dir+'gene_data.txt','w')
#output_log_handle = open(output_dir+'snp_parsing_log.tsv')

###log stuff
h1_snp_count = 0
h2_snp_count = 0
h3flag = False #boolian to use or not use third hypothesis
h3_snp_count = 0
n_a_count = 0
total_snp_count = 0
total_genes = 0

###Initiation sets
counts = [0,0,0,0,0,0] #Keeps track of the number of snps for each genome
#key_list = [] #allows reordering of dictionary (no longer used)
current_loci = '' #initiate gene parsing check
genome_list_for_highlighting = []
highlight_add_list = ['','','','','',''] # gather highlighting data in a bad way
highlight_add_list2 = ['','','','','','']

for line in input_file_handle:
    if line.startswith('LocusID'): #avoid header
        continue
    elements = line.split(',')
    if elements[31] == 'n/a': #avoid snps not in genes
        n_a_count += 1
        continue
    total_snp_count += 1
    # gather snp_data
    snps = elements[1:8]
    snp_location = elements[21]
    #other_stuff
    start, stop = elements[31].split('..'); start = int(start); stop = int(stop)
    loci = elements[32]; p = elements[27]; product = elements[33]
    chromosome = elements[0].split('.')[0]


    if current_loci == loci: #Add a count to each genome if it matches or doesnt match reference
        for s in range(6):# cycle through each snp
            if s != 0: # avoid checking reference
                if snps[s] == snps[0]: #add or subtract each one
                    counts[s]+=1
                else:
                    counts[s]-=1
    else:
        total_genes += 1
        #Write Gene Data
        for g in range(6):
            g += 1
            with open(output_dir+'genome'+str(g)+'_chr_'+chromosome[-1]+'.txt','a') as out:
                out.writelines(' '.join([chromosome,str(start),str(stop),str(counts[g-1]),'\n']))
                #print (g,chromosome, start, stop, counts[g-1])


        #Reset gene data
        current_loci = loci
        counts = [0,0,0,0,0,0]

        #grab initial gene data
        for s in range(6):# cycle through each snp
            if s != 0: # avoid checking reference
                if snps[s] == snps[0]: #add or subtract each one
                    counts[s-1]+=1
                else:
                    counts[s-1]-=1

    #####HIGHLIGHTING STUFF

    color = '' # color for highlighting, changes for each snp
    genome_list_for_highlighting = []

    #Hypotheis 1:XXOX*XO
    if p[0] == p[1] and p[1] == p[3]and p[1] == p[5] and p[1] != p[6] and p[1] != p[2]:# and p[1] != p[4]:
        genome_list_for_highlighting = [1,3,5]
        color = 'color=dred'
        h1_snp_count += 1
        if p[1] == p[4]:
            genome_list_for_highlighting.append(4)

    #Hypothesis 2: XXOX*XX
    elif p[0] == p[1] and p[1] == p[3] and p[1] == p[5] and p[1] == p[6] and p[1] != p[2]:# and p[1] != p[4]:  ### Light
 Orange
        genome_list_for_highlighting =[1,3,5,6] #ST36, 631(MAVP-Q), and 636
        color = 'color=yellow'
        h2_snp_count += 1
        if p[1] == p[4]:
            genome_list_for_highlighting.append(4)

    #Hypothesis 3: XXOOOX # only if you want
    elif p[1] == p[0] and p[1] == p[5] and p[1] == p[6] and p[1] != p[2] and p[1] != p[3] and p[1] != p[4]:
        if h3flag:
            genome_list_for_highlighting = [1,5,6]
            color = 'color=lorange'
        h3_snp_count += 1

## gather HIGHLIGHT OUTPUT for each genome
    for g in genome_list_for_highlighting:
        highlight_add = ' '.join([chromosome,snp_location,snp_location,"1",color,"\n"])
        if chromosome[-1] == '1':
            highlight_add_list[g-1]+=(highlight_add)
        else:
            highlight_add_list2[g-1]+=(highlight_add)

#### FOR GENE IDENTIFICATION
    if color:
        gene_info_handle.writelines(product + '::' +  loci + '::' + chromosome + '::' +  str(start) + '::' +  str(stop)
+ '::' +  color + '::' +  elements[34].rstrip()+'\n')


##Write highlighting info to output
for g in range(6):
    g += 1
    with open(output_dir+'genome'+str(g)+'_chr_'+'1.txt','a') as out:
        out.writelines(highlight_add_list[g-1])
for g in range(6):
    g += 1
    with open(output_dir+'genome'+str(g)+'_chr_'+'2.txt','a') as out:
        out.writelines(highlight_add_list2[g-1])






print ('Total number of snps:', total_snp_count)
print ('Total number of genes:', total_genes)
print ('Total snps for h1 (red):', h1_snp_count)
print ('Total snps for h2 (yellow):', h2_snp_count)
print ('Total snps for h3 (not shown):', h3_snp_count)

print ('Total snps with no gene information', n_a_count)
