#1)Viene lanciato un BLAST per individuare la regione genica in un set di genomi, e includendo oltre alla variante genica matchata con la query, una serie di basi a monte e a valle della regione. 

import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

tabella_blast=sys.argv[1].strip() #output di blast
genome=sys.argv[2].strip()  #genomi
chomp_size=int(sys.argv[3].strip())  #numero di tagli di basi

tabella = pd.read_csv(tabella_blast, encoding='utf7')

name = tabella["genome"]
start = tabella["start"]
end = tabella["end"]

outname="%s_chomped%s" %(genome, chomp_size)
outF=open(outname, "w")


inF=open(genome, "r")


for record in SeqIO.parse(genome, "fasta"):
	#inF=open(genome, "r")	
	ID=str(record.id)	
	for i in range(len(name)):	
		if ID == name[i]:   #riconosce numero ID nel confronto del for successivo		
			START=int(str(start[i]))			
			END=int(str(end[i]))			
										
			PRE_gdhB=START - chomp_size 
			POST_gdhB=END + chomp_size
			monte=str(record.seq[int(PRE_gdhB):START])
			valle=str(record.seq[END:int(POST_gdhB)])			
							
			if PRE_gdhB<0:
				monte=str(record.seq[:START])
			if len(valle) < chomp_size:
				valle=str(record.seq[END:])
			#print monte + '\n' + valle		
			if START > END:
				#START, END=END, START
				#5789,4347 
				
				PRE_gdhB=START + chomp_size 
				POST_gdhB=END - chomp_size
											
				monte=str(record.seq[START:int(PRE_gdhB)].reverse_complement())					
				valle=str(record.seq[int(POST_gdhB):END].reverse_complement())
				#print monte + '\n' + valle				
				if len(monte) < chomp_size:
					monte=str(record.seq[START:].reverse_complement())
				if POST_gdhB<0:
					valle=str(record.seq[:END].reverse_complement())				
											
			regione_monte=SeqRecord(Seq(monte), id=ID+"_PREgdhB", name="", description="")
			regione_valle=SeqRecord(Seq(valle), id=ID+"_POSTgdhB", name="", description="")
					
			regione_monte_valle=regione_monte, regione_valle			
						
			SeqIO.write(regione_monte_valle, outF, "fasta")
			

inF.close()
outF.close()


#cmd="cat *_chomped5000 >multi_chomped.fna | prodigal -i multichomped.fna -a GenQuery.faa -d gdhB1.gene | makeblastdb -in GenQuery.faa -dbtype prot -o GenQuery_db"
#os.system(cmd)





          





#cmd="for i in *.fasta;do makeblastdb -in $i -dbtype nucl; done"  
#os.system(cmd)
#cmd="for i in *.fasta;do blastn -db $i -query gdhB_partial -outfmt '6 sseqid sstart send perc_identity 95' -out $i; done"
#os.system(cmd)


