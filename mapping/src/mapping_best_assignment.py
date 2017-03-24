import argparse
import sys
import math
import re
from operator import itemgetter

parser = argparse.ArgumentParser(description='Find best barcode-locus assignment.')
parser.add_argument('-q', required=False, default=20, type=int, help="minimum mapping score (default 20)")
parser.add_argument('--min-reads', required=False,type=int,default=10,help="minimum # of reads to validate an assignment (default 10)")
parser.add_argument('--min-score', required=False,type=int,default=10,help="minimum assignment score (default 10)")
parser.add_argument('mapping_sam', help="SAM alignment file of integration loci with barcode in sequence name")
params = parser.parse_args()

bwa_q_thr = params.q
min_reads = params.min_reads
min_score = params.min_score
binsize = 100
strand = ['-','+']
maps = {}
bcdset = set()


with open(params.mapping_sam) as f:
   for line in f:
      if line[0] == '@': continue
      m = line.split('\t')
      # Exclude self-ligated reads.
#      if m[2] == '*' or m[2] == 'HIV' or int(m[4]) < bwa_q_thr: continue
      if m[2] == '*' or int(m[4]) < bwa_q_thr: continue
      # Parse SAM fields
      barcode = m[0]
      rchr = m[2]
      rloc = int(m[3])
      flags = int(m[1])
      dirflag = (flags & 16) > 0
      cigar = m[5]

      # Don't allow soft clipping at the LTR end.
      asymb = re.split('[0-9]+',cigar)[1:]
      avals = re.split('M|S|I|H|D',cigar)[:-1]
      if (dirflag and asymb[-1] == 'S') or (not dirflag and asymb[0] == 'S'):
         continue

      # Correct forward strand locus offset.
      if dirflag:
         for i,symb in enumerate(asymb):
            if symb == 'M' or symb == 'D':
               rloc += int(avals[i])
         
      # Add barcodes to loci dict
      keypos = rloc/binsize
      insert = 1
      for k in range(keypos-1,keypos+2):
         key = strand[dirflag] + rchr + str(k)
         if maps.has_key(key):
            insert = 0
            bcds = maps[key][1]
            if bcds.has_key(barcode):
               bcds[barcode] += 1
            else:
               bcds[barcode] = 1
            break
      if insert:
         key = strand[dirflag] + rchr + str(keypos)
         maps[key] = [rchr+":"+str(rloc)+":"+strand[dirflag],{barcode:1}]
         
      bcdset.add(barcode)

#DEBUG
"""      
print "loci-barcode table:"
for loc in maps:
   sys.stdout.write(maps[loc][0])
   bcds = maps[loc][1]
   for bcd in bcds:
      sys.stdout.write('\t('+str(bcds[bcd])+')'+bcd)
   sys.stdout.write('\n')
"""

# Convert bcd set in list
bcdlist = list(bcdset)
bcdsum  = [0]*len(bcdlist)
del bcdset
# Sort barcodes
bcdlist.sort()
# Create loc list
loclist = maps.keys()
locsum  = [0]*len(loclist)
# Build edge matrix
bcdcnt = len(bcdlist)
loccnt = len(maps)
#edges = [[0 for x in range(bcdcnt)] for x in range(loccnt)]
edges = []

# Make edge list
for l,loc in enumerate(loclist):
   bcds = maps[loc][1]
   for bcd in bcds:
      # Add read count to locus
      locsum[l] += bcds[bcd]
      # Add read count to barcode
      i = bcdlist.index(bcd)
      bcdsum[i] += bcds[bcd]
      # Assign read count to matrix
      edges.append([i,l,bcds[bcd]])
   del maps[loc][1]

# Sort edges.
edges = sorted(edges,key=itemgetter(2),reverse=True)

for edge in edges:
   bcd = edge[0]
   loc = edge[1]
   cnt = edge[2]
   if locsum[loc] < 0 or bcdsum[bcd] < 0: continue
   score = max((1-cnt*1.0/bcdsum[bcd]),(1-cnt*1.0/locsum[loc]))
   if score != 0:
      score = int(round(max(0,-10*math.log10(score))))
   else:
      score = 100
   if cnt >= min_reads and score >= min_score:
      sys.stdout.write(bcdlist[bcd]+'\t'+maps[loclist[loc]][0].replace(':','\t')+'\t'+str(cnt)+'\t'+str(score)+'\n')
      bcdsum[bcd] = -1
      locsum[loc] = -1

for loc,cnt in enumerate(locsum):
   if cnt >= min_reads:
      locus = maps[loclist[loc]][0].split(':')
      if locus[0] == 'HIV':
         continue
      sys.stdout.write('NA\t'+maps[loclist[loc]][0].replace(':','\t')+'\t'+str(cnt)+'\tNA\n')
