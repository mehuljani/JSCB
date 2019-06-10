
import sys 
import os
import random
from Bio import SeqIO

if len(sys.argv) == 1:
	print "\n\nError: Please provide input file in genbank format!!\n\n"
	sys.exit()
f = open ("JSCBinput", "w")
start_col, end_col, rna_start, rna_end, strand_col=[],[],[],[],[]
actstart = []
actend = []
nucs = ['A','T','C','G']
for gb_record in SeqIO.parse(sys.argv[1], "genbank"):
	seq=gb_record.seq
	seq=str(seq.upper())
	finseq=""
	for i in seq:
		if i in nucs:
			finseq+=i
		else:
			finseq+=random.choice(nucs)

	f.write(finseq)
	f.close()
	
	for feature in gb_record.features:
		start=feature.location.nofuzzy_start+1
		end=feature.location.nofuzzy_end
		if feature.strand == 1:
			strand ="+"
		else:
			strand="-"
		if feature.type == "CDS":
			start_col.append(start)
			end_col.append(end)
			strand_col.append(strand)
		if feature.type == "rRNA":
			rna_start.append(start)
			rna_end.append(end)
	
	f1 = open ("JSCB_coord", "w")		
	for x, y, z in zip(strand_col,start_col,end_col):
		if (z-y+1)%3==0:
			f1.write(str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
			actstart.append(y)
			actend.append(z)
	f1.close()			

os.system("./jscb")

gifile = open("JSCB_output.gi", "r")
hashclussize = {}
for line in gifile:
	try:
		hashclussize[line.split()[1]]+=1
	except:
		hashclussize[line.split()[1]]=1

natclus = max(hashclussize, key=hashclussize.get)

prv = -1
gmode = 0
gicounter=0
gifile = open("JSCB_output.gi", "r")
outfile = open("JSCB_output.tsv", "w")

for line in gifile:
	if int(line.split()[1])!=int(natclus):
		if int(line.split()[0]) != prv+1:
			startgeneno = int(line.split()[0])
			startcoord = actstart[int(line.split()[0])-1]
			gmode = 1
		prv = int(line.split()[0])
	elif gmode == 1:
		endgeneno = int(line.split()[0])
		gilength = endgeneno-startgeneno
		endcoord = actend[int(line.split()[0])-2]
		if gilength >= 8:
			rnais = 0
			for x, y in zip(rna_start, rna_end):
				if startcoord < x and startcoord < y and endcoord > x and endcoord > y:
					rnais = 1
				
			if rnais == 0:
				gicounter+=1
				if gicounter == 1:
					outfile.write("Id\tStart\tEnd\n")
				outfile.write("GI-"+str(gicounter)+"\t"+str(startcoord)+"\t"+str(actend[int(line.split()[0])-2])+"\n")					
				
		gmode = 0
if gicounter == 0:
	outfile.write("No Genomic Islands Identified\n")
outfile.close()
os.system("rm JSCBinput")
os.system("rm JSCB_coord")
