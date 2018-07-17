# -*- coding: utf-8 -*-
import os
import sys
import numpy
import string
import commands
from MyClasses import *
from optparse import OptionParser


parser = OptionParser(usage = "%prog -f FILE -p FILE -d FILE")
parser.add_option("-d", "--dot", dest = "dotf", metavar = "FILE" ,help = "dot file")
parser.add_option("-p", "--pdb", dest = "pdbf", metavar = "FILE",help = "PDB attributes")
parser.add_option("-f", "--fas", dest = "fasf",  metavar = "FILE" ,help = "fasta file")
parser.add_option("-P", "--path", dest = "path",  metavar = "PATH" ,help = "path to programs")
(args, options) = parser.parse_args()
dotfile = args.dotf
#path for all Programs
#cheminPro = "./Programs"
dotf = dotfile.split('/')[1].split('-col.')[0]
fastafile = args.fasf
def main(dotfile = dotf, attrfile = args.pdbf, fastaFile = args.fasf, cheminPro = args.path):
        chemin = fastafile.replace('.fasta', "_paloma")
	
	sequencesdep = getSequencesOf(fastaFile)
	if os.path.exists(fastafile.replace(".fasta","_r3D")):
		commands.getoutput('rm -r '+fastafile.replace(".fasta","_r3D")+"/"+dotfile)
	commands.getoutput('mkdir '+fastafile.replace(".fasta","_r3D")+"/"+dotfile)
	fichierr = open( chemin+"/"+dotfile+"-col.dot", "r")
	fichierw = open( dotfile+'.xml', "w")
	fichierw.write('''<?xml version=\"1.0\"?><resume>''')

	sequences, sequences2, dico_plma, dico_plma2 = [], [], dict(), dict()
	dico_plma_seq, loc, sub, ex = dict(), 0, 0, 0
	########################################################################################################
	for line in fichierr:
		if 'subgraph cluster_1\n' in line: loc = 1
		elif 'subgraph cluster_2\n' in line:loc= 2
		elif 'subgraph cluster_3\n' in line:loc= 3
		elif '->' in line: loc=4

		if loc==1: #pour les noms de sequences
			if line.strip()[0]=='\"':
				a = line.split(':')[1]
				sequences2.append(line.split(':')[1].split('\"')[0])
				sequences.append([])
	
		if loc==3: #pour les plmas
			if 'subgraph cluster' in line:
				sub=line[line.find('_')+1:-1]
				dico_plma['Plma'+sub] = []
				dico_plma2['Plma'+sub] = []
				
			if line.strip()[0]=='\"':
				ip = line.split('\"')
				
				seq = int(ip[1].split(',')[0][1:])
				pos = int(ip[1].split(',')[1][1:])
				length = int(ip[1].split(',')[2][1:-1])
				plma = ip[3]
				sequences[seq-1].append((pos,'Plma'+sub))
				dico_plma_seq['Plma'+sub+'_'+str(seq)] = (pos, pos+length)

				for j in range(length):
					if len(dico_plma['Plma'+sub]) != length: 
						dico_plma['Plma'+sub].append([(seq,pos)])
						dico_plma2['Plma'+sub].append([plma[j]])
					else:
						dico_plma['Plma' + sub][j].append((seq, pos))
						if plma[j] not in dico_plma2['Plma'+sub][j]:
                                                        dico_plma2['Plma'+sub][j].append(plma[j])		
		if loc==4: #pour les fleches au niveau des couleurs
			if line.strip()[0]=='\"':
				ip = line.split('\"')
				seq = int(ip[1].split(',')[0][1:])
				pos = int(ip[1].split(',')[1][1:])
				lg = int(ip[1].split(',')[2][1:-1])
				pos2 = int(ip[3].split(',')[1][1:])-1
				if 'style' in line: 
					ex = ex+1
					sequences[seq-1].append((pos+1,'PlmaEx'+str(ex)))
					l = []
					for p in sequencesdep[seq-1][pos-1+lg:pos2]:
						l.append([p])
					dico_plma2['PlmaEx'+str(ex)] = l
					
					
				elif '..' in line: sequences[seq-1].append((pos+1,'Gap'))
        #######################################################################################################
	newseq = []
	for iSeq in sequences:
		part = []
		for j in sorted(iSeq, key = lambda x: x[0]):
			part.append(j[1])
		newseq.append(part)
	
	diff = fileconcept(newseq, sequences2, dico_plma2, fastafile.replace(".fasta","_r3D")+"/"+dotfile)

	pdb = getpdbinfo(attrfile)

	fileconc = open(os.getcwd()+"/concepts","w")
	conceptsint = []

        #######################################################################################################
	for seq,p in pdb.iteritems():
		fileconc.write(p[0]+"_2.pdb\n")
		for plma in newseq[int(seq)-1]:
                        
			if plma+"_"+str(int(seq)) in dico_plma_seq: 
				fichiertestpdb = open(os.getcwd() +"/"+fastafile.replace(".fasta","_r3D")+"/"+dotfile +"/"+plma, "w")
				fichiertestpdb.write(">"+sequences2[int(seq)-1]+"\n")
				fichiertestpdb.write(sequencesdep[int(seq)-1]+"\n")
				
				for plma2 in diff[0][plma]: 
					fichiertestpdb.write(str(dico_plma_seq[plma2 + "_" + str(int(seq))][0])
                                             + "-"+ str(dico_plma_seq[plma2+"_" + str(int(seq))][1])+"-blue\n")
					
				for plma2 in diff[1][plma]: 
					fichiertestpdb.write(str(dico_plma_seq[plma2 + "_" + str(int(seq))][0])
                                             + "-" + str(dico_plma_seq[plma2+"_"+str(int(seq))][1]) + "-red\n")
					
				fichiertestpdb.write( str( dico_plma_seq[plma + "_" + str(int(seq) )][0]) + "-"
                                              + str(dico_plma_seq[plma + "_" + str(int(seq))][1]) + "-green\n")
				fichiertestpdb.close()
				print commands.getoutput("python "+cheminPro+"/alignpdbfas.py -f " +fastafile.replace(".fasta","_r3D")+"/"+dotfile
				+"/"+plma+" -p ./PDB/"+p[0]+".pdb -c "+p[1]+" -o test.html -O ./PDB/"+p[0]+".pdb")
				
				colorpdb(os.getcwd()+'/PDB/'+p[0]+'.pdb',p[1],os.getcwd()+'/'+fastafile.replace(".fasta","_r3D")+"/"+dotfile
                            +"/"+plma,os.getcwd()+"/"+fastafile.replace(".fasta","_r3D")+"/"+dotfile +"/"+ plma+".py", "./PDB/"+p[0]+".pdb")
				 
	
	fileconc.write("\n\n\n Interessting blocks")
	conceptsint2 = sorted(conceptsint, key = lambda x: x[2])
	#######################################################################################################
	for z in enumerate(conceptsint2): 
		fileconc.write("\n")
		for y in z[1]: 
			if z[0]%2 == 0:
				fileconc.write(str(y)+" ")
		
	fileconc.close()
	num = 1	
	fichierw.write('\t<sequences>\n')
	######################################################################################################
	for i,j in enumerate(sequences):
		fichierw.write('\t\t<sequence tag=\"'+str(num)+": "+sequences2[i]+'\">\n')
		num = num + 1
		for k in sorted(j,key=lambda x: x[0]): 
			fichierw.write('\t\t\t<part data=\"'+k[1]+'\">\n')
			fichierw.write('\t\t\t</part>\n')
		fichierw.write('\t\t</sequence>\n')
	fichierw.write('\t</sequences>\n')

	fichierw.write('\t<blocks>\n')
	######################################################################################################
	for i,j in sorted(dico_plma2.iteritems()):
		fichierw.write('\t\t<block tag=\"'+i+'\">\n')
		for k in j:
			fichierw.write('\t\t\t<part>\n')
			for z in k:
				fichierw.write('\t\t\t\t<aa data=\"'+z+'\">\n')
				fichierw.write('\t\t\t\t</aa>\n')

			fichierw.write('\t\t\t</part>\n')
		fichierw.write('\t\t</block>\n')
	fichierw.write('\t</blocks>\n')
	fichierw.write('</resume>\n')
	makecliquabledot(fastafile.replace(".fasta","_r3D"),chemin+"/"+dotfile+"-col.dot") 


	###### recodage en gardant la longueur
	newsequences = []
	seq = []
	######################################################################################################
	for s1 in newseq: 
		for plma in s1:
			#print plma
			if plma[0] == 'G': seq.append(plma)
			else: 
				for i,pos in enumerate(dico_plma2[plma]): 
					seq.append(plma+'_'+str(i))
				
		newsequences.append(seq)
		seq = []
					
	newseq = newsequences
	dico_plma3 = dict()
	
	for plma,l in dico_plma2.iteritems():
		for i,pos in enumerate(l):
			dico_plma3[plma+'_'+str(i)] = [pos]
		
	dico_plma2 = dico_plma3		
	rec = (sequencesdep,newseq,dico_plma2)
	

###########################################################################################
def colorpdb(pdbfile,chain,blocs,outfile, pdbId):
        """
        This function makes it possible to write a file containing commands which are comprehensible
        and interpretable by the pymol software.
        When clicking on a given PLMA block, the pymol software is called. It loads and executes
        the commands contained in the file.
        """
	fichierw = open(outfile, "w")
	fichierw.write("from pymol import *\n")
	fichierw.write("cmd.load('"+ pdbfile +"')\n")
	fichierw.write("cmd.hide('everything')\n")
	fichierw.write("cmd.show('cartoon')\n")
	fichierw.write("cmd.color('gray')\n") 
	fichierw.write("cmd.bg_color('white')\n")

	fichier = open(blocs, "r")
	for i in fichier: 
		if i.strip()[0] in ['0','1','2','3','4','5','6','7','8','9']:
			a = i.strip().split('-');
			fichierw.write("cmd.color('"+a[2]+"', 'resi "+a[0]+"-"+a[1]+" and chain "+chain+"')\n")
	#zoom for displaying
	fichierw.write("cmd.zoom('chain "+ chain +"',0)\n") 	
			
        #The residues of the active site of the protein are recovered
	residues = getActiveSiteOf(pdbId)
	for line in residues:
                if line[1].strip()!='':
                        fichierw.write("cmd.show('sticks','resi "+line[1]+"')\n")
                        fichierw.write("cmd.color('orange','resi "+line[1]+"')\n")
        fichierw.write("cmd.ray\n")        
	fichier.close()
	fichierw.close()
############################################################################################################
def getActiveSiteOf(struct):
        """
        This function opens a pdb file and retrieves the catalytic region of the active site
        """
        aas=['ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL','XXX']
        commands.getoutput('grep ^SITE '+struct+ '> ./PDB/test.txt')
        with open ('./PDB/test.txt', 'r') as fd:
                struct_sites = fd.readlines()
                sites, i = [], 0
                while i < len(struct_sites):
                        line = struct_sites[i]
                        if line[8:10].strip() == str(1):
                                j=23
                                while j+11 < len(line):
                                        if line[j-5:j-2] in aas:
                                                chain_pos = []
                                                chain_pos.append(line[j-1:j])
                                                chain_pos.append(line[j+1:j+4])
                                                sites.append(chain_pos)
                                        j+=11
                                while i+1 < len(struct_sites) and struct_sites[i+1][8:10].strip() != str(1):
                                        line = struct_sites[i+1]
                                        j=23
                                        while j+11 < len(line):
                                                if line[j-5:j-2] in aas:
                                                        chain_pos = []
                                                        chain_pos.append(line[j-1:j])
                                                        chain_pos.append(line[j+1:j+4])
                                                        sites.append(chain_pos)
                                                j+=11
                                        i+=1
                                i+=1
        commands.getoutput('rm ./PDB/test.txt')
        return sites
#####################################################################################
def testprox(pdb,deb1,deb2,fin1,fin2):
	structure = myPDBParser(pdb,'firstChainOnly')
	try:
		rect = structure.getIndexFromResIDRect(((deb1,'A'),(fin1,'A'),(deb2,'A'),(fin2,'A')))
		dm = structure.getPortionOfDM(rect)
	except IndexError:
		return (False,0)
	except TypeError: 
		return (False,0)
	else:
		dist = numpy.amin(dm)
		if dist <= float(7):
			return (True,dist)
		else: return (False,0)

###############################################################################################################
def makecliquabledot(directory,dot): 
	fichierw = open(dot+"bis","w")
	fichier = open(dot, "r")
	node = False
	sub = 0
	for i in fichier:
                if "subgraph" in i.strip():
                        try:
                                x1 = int(i.strip().split('_')[-1])
                        except:
                                pass
                        plma = 'Plma'+i.strip().split('_')[-1]
		dotok = dot.split('/')[1]
		if "subgraph" in i.strip() and isinstance(x1, int) and x1 >2 and plma in os.listdir(os.getcwd() +"/"+directory+"/"+dotok.split('-col.')[0]):
			node = True
			sub = i.strip().split('_')[-1]
		elif node and "node" in i:
			node = False
			fichierw.write("target=\"file:////usr/bin/pymol\";\n")
			fichierw.write("URL=\"file://"+os.getcwd()+"/"+fastafile.replace(".fasta","_r3D")+"/"+dotok.split('-col.')[0]+"/Plma"+sub+".py\";\n")	
		fichierw.write(i)

	for i in os.listdir(os.getcwd()):
		if 'tmp' in i or '.xml' in i or i=='test.html': 
			commands.getoutput("rm -f "+i)

	fichier.close()
	fichierw.close()
	commands.getoutput("dot -Tps2 Gsize=\"100,100\" "+dot+"bis -o "+dot+"bis.ps")
	commands.getoutput("ps2pdf "+dot+"bis.ps "+dot+"3D.pdf")
	#commands.getoutput("acroread "+dot+"3D.pdf")
###############################################################################################################	
def getSequencesOf(filename):
	sequences, seq = [], ''
	fichier = open(filename, "r")
	for line in fichier: 
		if line[0] != '>': seq += line[0:-1]
		else:
                        sequences.append(string.upper(seq))
			seq = ''
	sequences.append(string.upper(seq))	
	fichier.close()
	return sequences[1:]
##############################################################################################################
def fileconcept(newseq, seq, dicoplma, dotfile):
	dicoId, dico_seq, Id = {}, {}, 0
	
	for i,j in dicoplma.iteritems(): 
		if 'Ex' not in i: 
			dico_seq[i] = []
			dicoId[i] = Id
			Id = Id + 1
			
	for i,j in enumerate(newseq): 
		for plma in j:
			if 'Ex' not in plma and 'Gap' not in plma: 	
				dico_seq[plma].append(i)
				
	return differentiel(dico_seq,newseq)
#############################################################################################################
def differentiel(dico, seq):
	dicodiff, dicoglob = {}, {}
	for i,j in dico.iteritems():
		dicodiff[i], lp = [], []
		for sequences in j: 	
			lp.append(set(seq[sequences]))
			
		globale = list(set.intersection(*lp))
		if 'Gap' in globale: globale.remove('Gap')
		dicoglob[i] = globale[:]
		
		for k in globale: 
			xx =  set(dico[k])
		 	if set(j) == xx:
				 dicodiff[i].append(k)
				 dicoglob[i].remove(k)

		if i in dicodiff[i]: dicodiff[i].remove(i)
		
	return (dicoglob,dicodiff)

###########################################################################################################
def getpdbinfo(filepdb): 
	from xml.dom import minidom
	filein = filepdb
	dico = {}
	dom = minidom.parse(filein)
	sequences = dom.getElementsByTagName('sequences')
	partition = sequences[0].getElementsByTagName('seq')
	for part in partition:
		dico[part.getAttribute("id")]=(part.getAttribute("pdb"),part.getAttribute("chain"))

	#for i,j in dico.iteritems(): 
		#print i,j[0],j[1]	

	return dico 

##########################################################################################################
if __name__ == '__main__':
	main()
	

