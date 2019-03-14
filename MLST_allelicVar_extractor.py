# Libraries
import os
import sys
from Bio import SeqIO
import pandas as pd




# main
### fetch input

gFile=sys.argv[1].strip()

sFile=sys.argv[2].strip()

dFile=sys.argv[3].strip()

####Error file
outname="%s_error.log" %(gFile)
outF=open(outname,"w")



### blast sequences on genome

# qseqid sseqid pident qlen length | qstart qend slen sstart send

cmd="makeblastdb -in %s -dbtype nucl" %(gFile)
os.system(cmd)

cmd="blastn -query %s -db %s -outfmt '6 qseqid sseqid pident qlen length qstart qend slen sstart send' -out %s_MLSTout -perc_identity 98" %(sFile,gFile,gFile)
os.system(cmd)

### parse output

FF="%s_MLSTout" %(gFile)
inF=open(FF,"r")

hash={}
for i in inF.xreadlines():
	i=i.strip().split()
	if (float(i[2])==100) and (i[3]==i[4]):
		gene=i[0].strip().split("_")[:-1]
		gene="_".join(gene)
		var=i[0].strip().split("_")[-1]
		if gene in hash.keys():
			outF.write("more than one allele for gene %s in genome %s\n" %(gene,gFile))
		else:
			hash[gene]=var


inF.close()

if len(hash.keys())!=7:
	outF.write("error in allele extraction for file %s\n" %(gFile))
	for y in hash.keys():
		outF.write("%s\t%s\n" %(y,hash[y]))

else:

	df=pd.read_csv(dFile, sep="\t")

	dm=df.copy()
	for x in hash.keys():
		#print x
		#print hash[x]
		dm=dm.loc[dm[x] == int(hash[x])]
	
	if len(dm)!=1:
		outF.write("error in MLST definition for file %s\n" %(gFile))
		for y in hash.keys():
			outF.write("%s\t%s\n" %(y,hash[y]))
	
	else:
		csv_name="%s_result.tsv" %(gFile)
		dd=dm.iloc[:,0:8]
		dd = dd.assign(genome=gFile)
		dd.to_csv(csv_name,sep="\t",header=False, index=False)

outF.close()

