import sys
import seeq

seqno = 1
bcd_lookup = [""]*100000000

with open(sys.argv[1]) as fs, open(sys.argv[2]) as fq:
   for line in fs:
      tok = line.split("\t")
      bcd = tok[0]
      ids = tok[2].split(",")
      if not bcd:
         bcd = "NA"
      for id in ids:
         bcd_lookup[int(id)] = bcd

   for lineno,line in enumerate(fq):
      if lineno%4 == 0:
         sys.stdout.write("@"+bcd_lookup[seqno]+"\n")
         seqno += 1
      else:
         sys.stdout.write(line)
