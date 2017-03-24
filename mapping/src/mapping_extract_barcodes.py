import argparse
import seeq
import gzip

parser = argparse.ArgumentParser(description='Extract barcode and genome sequence from iPCR reads.')
parser.add_argument('-b',required=True,help="barcode flanking sequence (at 3')")
parser.add_argument('-l',required=True,help="construct sequence present in genome (reverse) reads")
parser.add_argument('-r',required=False,default='',help="restriction enzyme flanking sequence (at construct 3')")
parser.add_argument('--bdist', required=False,type=int,default=-1,help="matching distance for -b. (Default max(1,0.2*length))")
parser.add_argument('--ldist', required=False,type=int,default=-1,help="matching distance for -l. (Default max(1,0.1*length))")
parser.add_argument('--rdist', required=False,type=int,default=0,help="matching distance for -r. (Default 0)")
parser.add_argument('--barcodes', default='ipcr.bcd', help='barcode output file')
parser.add_argument('--reads', default='reads.fastq', help='genome reads output file')
parser.add_argument('ipcr_forward', help="iPCR forward reads (gzipped+fastq format)")
parser.add_argument('ipcr_reverse', help="iPCR reverse reads (gzipped+fastq format)")
params = parser.parse_args()

fbrcd = open(params.barcodes, 'w')
fseqs = open(params.reads, 'w')

if params.bdist < 0: bdist = int(max(1,round(0.2*len(params.b))))
else:                bdist = params.bdist
if params.ldist < 0: ldist = int(max(1,round(0.1*len(params.l))))
else:                ldist = params.rdist

# BHIVE seqs:
# T7 promoter (-b) TATAGTGAGTCGTA
# LTR sequence (-l) AGCCCTTCCA
# HIVRE sequence (-r) CGCTTTTAA
T7 = seeq.compile(params.b,bdist)
LTR = seeq.compile(params.l,ldist)
HIVRE = None
if params.r:
   HIVRE = seeq.compile(params.r,params.rdist)
fqline = [""]*4

with gzip.open(params.ipcr_forward) as r1, gzip.open(params.ipcr_reverse) as r2:
   for lineno,line in enumerate(r1):
      if lineno % 4 != 1: continue
      for i in range(0,4):
         fqline[i] = r2.readline()
      # Match 5' LTR on reverse read (required, there must be 20 nt).
      l = LTR.matchBest(fqline[1])
      if not l:
         continue

      # Match HIV BplI site (not required)
      m = HIVRE.match(fqline[1]) if HIVRE else None
      beg = l.matchlist[0][1]
      end = len(fqline[1])
      if m:
         end = m.matchlist[0][0]
         
      # Update fastq sequence
      if end - beg < 1:
         continue
      fqline[1] = fqline[1][beg:end].rstrip()+"\n"
      fqline[3] = fqline[3][beg:end].rstrip()+"\n"

         
      # Match T7 promoter (not required)
      brcdline = "\n"
      m = T7.matchBest(line)
      if m:
         brcd = m.tokenize()[0]
         brcdline = brcd + '\n'

      # Write barcode
      fbrcd.write(brcdline)
      
      # Write fastq read
      for i in range(0,4):
         fseqs.write(fqline[i])
      
r1.close()
r2.close()
fbrcd.close()
fseqs.close()
