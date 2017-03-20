#!/usr/bin/python3
import sys

#blastp -query ../10329.faa -db cog_db -evalue 1e-10 -num_threads 20 -outfmt '6 sseqid qseqid slen qlen qcovs pident len
gth qstart qend sstart send evalue bitscore' -out 10329_proteins_vs_cog_1e-10.blast

blast = sys.argv[1]
best_hits = {}



for line in open(blast,'r').readlines():
    elements = line.rstrip().split('\t')
    #seq_id = elements[0]; qseq_id = elements[1];
    sseqid, qseqid, slen, qlen, qcovs, pident, length, qstart, qend, sstart, send, evalue, bitscore = elements[:]
    #print (sseqid, qseqid, length, slen, qlen, pident)
    if qseqid in best_hits.keys():
        if float(bitscore) > float(best_hits[qseqid][1]):
            best_hits[qseqid] = [sseqid,bitscore,qlen,slen,length,pident]
    else:
        best_hits[qseqid] = [sseqid,bitscore,qlen,slen,length,pident]


for q,elements in best_hits.items():
    print (q, ' '.join(elements))
