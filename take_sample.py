	# importing module
from pandas import *
import sys
import numpy as np
import csv


#1) chemin du fichier contenant les echantillons et type 
#2) chemin du fichier contenant l'ensemble des matrices
#3) -o output 

#exemple python3 take_sample.py /path/samplehuf.test0.tsv /path/mat3_mer.tsv /path/outdir

sample_path=sys.argv[1]
mat_path=sys.argv[2]
outfile=sys.argv[3]

def take_name(sample_path):
#Fonction qui prend les identifiant contenut au niveau des samples
	list=[]
	lst_pos=[]
	name=[]
	position=-1
	file=open(sample_path,"r")
	for line in file:
			
		for char in line:
			position=position+1
			if char == '\t':
				lst_pos.append(position)
				name.append(line[0:position])
			elif char == '\n':
				position=-1	
		
	#for i in lst_pos :
	#	list.append(echant[0:i])				
	file.close()
	return name
name=take_name(sample_path)

#récupere les colonnes dans la matrices total en comparant le header avec les echantillons choisit
def take_sample(mat_path,name,outfile):
	#csv.field_size_limit(2000000)
	file = open(mat_path)
	reader = csv.reader(file, delimiter = '\t')
	indice=[]
	fichier=open(outfile,"a")
	for line in reader:
	#	print(line[1])
		if line[0]=="":
			i=0	
			while i < len(line):
				j=0
				while j < len(name):
					if name[j]==line[i]:
						print(name[j],line[i])
						indice.append(i)
					j=j+1
				i=i+1	
			fichier.write("{0}\t".format(line[0]))
		else:	
			fichier.write("\n{0}\t".format(line[0]))			
		for i in indice:
			fichier.write("{0}\t".format(line[i]))	
	return indice
test2=take_sample(mat_path,name,outfile)
print(test2)	
	
