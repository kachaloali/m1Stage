#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on wed May 24 08:40:31 2017

@author: Kachalo Ali
"""

import os
import Bio
import sys
import math
import string
import commands
from Bio import pairwise2
from Bio.PDB import PDBList
from optparse import OptionParser
from collections import OrderedDict
from Bio.pairwise2 import format_alignment


parser = OptionParser(usage="%prog -f FILE")
parser.add_option("-f", "--fasta", dest = "fasta", metavar = "FILE",
help ="fasta files and two files serving as databases separated by commas")
parser.add_option("-p", "--path", dest = "path", metavar = "PATH",
help ="path to programs")
(args, options) = parser.parse_args()

#main
def main(fastaFile = args.fasta.split(','), pathToProg=args.path):
    #the identifiers of the proteins are recovered from the fasta file
    diIdSeqs = getProtIds(fastaFile[0])
    #fastaFile[1] & fastaFile[2] represent the 1st and the 2nd file serving as databases
    diProt = getpdbIds(diIdSeqs,fastaFile[0], fastaFile[1],fastaFile[2], pathToProg)
    #the attribute file is created in this step
    getOutPut(diProt, fastaFile[0])

def getOutPut(diProt, fastaFile):
    """
    Allows to set the pdb attribute file associated with the input fasta file
    """
    pdb_attr = open(fastaFile.replace('.fasta','')+'_pdb.attr', 'w')
    pdb_attr.write('<sequences>\n')
    numSeq = 0
    for idProt, associate_pdb in diProt.items():
        if associate_pdb != '':
            pdb_attr.write('\t<seq id="'+str(numSeq+1)+'" pdb="'+associate_pdb+'" chain="A" />\n')
            numSeq +=1
        else: numSeq +=1
            
    pdb_attr.write('</sequences>')
    pdb_attr.close()

    
def getpdbIds(dic,fastafile, fbase1, fbase2, pathToProg):
    #old files are deleting 
    if os.path.exists("./PDB/"):
        commands.getoutput('rm  ./PDB/*.pdb')
    #create a new folder that will contain the pdb files to download
    commands.getoutput('mkdir ./PDB/')
    #file 1
    with open (fbase1, 'r') as fd:
        database1 = fd.readlines()
    #file 2    
    with open (fbase2,'r') as fd:
        database2 = fd.readlines()
    #Output dictionary
    diProt = {}
        
    pdbl = PDBList()
    '''Selecting structures from PDB'''
    for protId in dic:
        #The 3D structures associated with the protId protein in the first file
        associate1 = list(filter(lambda line: protId in line, database1))
        #The 3D structures associated with the protId protein in the second file
        associate2 = list(filter(lambda line: protId in line, database2))
        associate_pdb = []
        
        if associate1 != []:
            #When associated 3D structures are found in file 1
            for line in associate1:
                pdbId = line.split('\t')[0]
                associate_pdb.append(pdbId)
            pdbId = getBestPdbId(fastafile, protId, associate_pdb, pathToProg)
            if commands.getoutput('wget -c http://www.rcsb.org/pdb/files/'+pdbId+'.pdb -O '+os.getcwd()+"/PDB/"+pdbId+'.pdb'):
                diProt[protId] = pdbId
            else: diProt[protId] = ""
            
        elif associate2 != []:
            #When associated 3D structures are found in file 2
            for line in associate2:
                pdbId = line.split('\t')[0]
                associate_pdb.append(pdbId)
            pdbId = getBestPdbId(fastafile, protId, associate_pdb)            
            if commands.getoutput('wget -c http://www.rcsb.org/pdb/files/'+pdbId+'.pdb -O '+os.getcwd()+"/PDB/"+pdbId+'.pdb'):
                diProt[protId] = pdbId
            else: diProt[protId] = ""
                
        else: diProt[protId] = ""
    #Everything is renaming from .cif to .pdb if it exists
    commands.getoutput('bash ./Programs/rename.sh cif pdb')
    if os.path.exists("pdb_tmp"):
        commands.getoutput('rm -r ./pdb_tmp') 
    #Return a dictionary whose keys are the identifiers of the proteins
    #and the aqssociated values are the pdb identifiers associated with its proteins
    return diProt
      
def getProtIds(fasta):
    """
    Given a fasta file, the identifiers of the proteins will be recover
    """
    ifile =  open(fasta, 'r')
    infile = ifile.readlines()
    diIdseqs = {}
    for line in infile:
        if line[0] == '>' :
            if '|' in line:
                
                iSeq = line.split('|')[1].strip()
                if '_' in iSeq.strip(' '):
                    idSeq = iSeq.split('_')[0]              
                    chainSeq = iSeq.split('_')[1]
                    
                elif ':' in iSeq:
                    idSeq = iSeq.split(':')[0]              
                    chainSeq = iSeq.split(':')[1]
                    
                else:
                    idSeq, chainSeq = iSeq, ''           
                diIdseqs[idSeq] = chainSeq
                
            else:
                idSeq = line[1:].strip()
                
    ifile.close()
    return diIdseqs
#========================================================================================
#========================================================================================
def pdb2fasta(chemin, pdbId, pathToProg):
    """Allows to convert a 3D structure into amino acid sequences"""
    print commands.getoutput("./"+ pathToProg +"/pdb2fasta.sh " + chemin + " " + pdbId)
    return open(chemin[:-8] + pdbId +'.fasta', 'r').read().split('\n')[1].strip(), pdbId

        
def makeFasta(idSeq, descs, sequences):
    """Builds a fasta file for a single sequence"""
    i= 0
    while(i < len(descs)):
        if(idSeq in descs[i]):
            prot = open(idSeq + '.fasta', 'w')
            prot.write('>' + descs[i] + '\n')
            prot.write(sequences[i])
            prot.close()
            return sequences[i], idSeq
        i += 1


def parse_fasta (lines):
    """Returns 2 lists: one is for the descriptions and other one for the sequences"""
    descs, seqs, data = [], [], ''
    for line in lines:
        if line.startswith('>'):
            if data:   # have collected a sequence, push to seqs
                seqs.append(data)
                data = ''
            descs.append(line[1:])  # Trim '>' from beginning
        else:
            data += line.rstrip('\r\n')
    # there will be yet one more to push when we run out
    seqs.append(data)
    return descs, seqs


def parse_blastp_And_Get_Evalue (ifile,pdbId, score, m, n):
    """Parse the output of blastp to recover the e-values"""
    lines, i = open(ifile, 'r').read().split('\n'), 0
    #retrieving of the e-value
    while (i <=len(lines)):
        if lines[i].strip().startswith(pdbId):
            var = lines[i]
            eValue = float(var.split()[2])
            i += len(lines)
        i +=1
    #parse the output of blastp to recover the values: lambda & k
    #after that we compute the e-value of the alignment to be recovered
    if i != len(lines):
	#Case where the sequences are very short in the alignment file
        #the e-value is not displayed, it will have to be computed
        i = 0
        while (i <=len(lines)):
            if "Lambda" in lines[i]:
                var = lines[i+1]
                #The e-value of the alignment is computed
                eValue = float(var.split()[1]) * m * n * math.exp(-float(var.split()[0]) * float(score))
                i += len(lines)
            i +=1
        
    return eValue
            

def getBestPdbId(ifile,sprotId,tab, pathToProg):
    """Selects the best pdb structures that must be associated"""
    dipdbAndEvalues = dict()
    commands.getoutput('mkdir ./pdb_tmp')
    lines = open(ifile, 'r').read().split('\n')
    descriptions, sequences = parse_fasta(lines)
    seq1 = makeFasta(sprotId, descriptions, sequences)
    commands.getoutput('mv '+sprotId+'.fasta ./pdb_tmp/')
     
    for pdbId in tab:
        if commands.getoutput('wget -c http://www.rcsb.org/pdb/files/'+pdbId+'.pdb -O ./pdb_tmp/' +pdbId+'.pdb'):
            seq2 = pdb2fasta("./pdb_tmp/"+ pdbId +".pdb", pdbId, pathToProg)
            print commands.getoutput("blastp -query ./pdb_tmp/" +seq1[1]+".fasta -subject ./pdb_tmp/"
                                                                    +seq2[1]+".fasta > ./pdb_tmp/blastp.txt")
            alignments = pairwise2.align.localxx(seq1[0], seq2[0])
            score = pairwise2.format_alignment(*alignments[0]).split('Score=')[1].strip()
            dipdbAndEvalues[pdbId] = parse_blastp_And_Get_Evalue('./pdb_tmp/blastp.txt',pdbId,score,len(seq1),len(seq2))
            if os.path.exists("/pdb_tmp/"+ pdbId +".pdb"):
                print commands.getoutput("rm ./pdb_tmp/"+ pdbId +".pdb")
            if os.path.exists("./pdb_tmp/"+ sprotId +".pdb"):
                print commands.getoutput("rm ./pdb_tmp/"+ sprotId +".pdb")
	    print commands.getoutput("rm ./pdb_tmp/"+ pdbId+".fasta")
        else: dipdbAndEvalues[pdbId] = ""
    if os.path.exists("./pdb_tmp/"+ sprotId+".fasta"):
        print commands.getoutput("rm ./pdb_tmp/"+ sprotId+".fasta")
    dipdbAndEvalues = OrderedDict(sorted(dipdbAndEvalues.items(), key = lambda t: t[1]))
    return dipdbAndEvalues.items()[0][0]
    

#main
if __name__ == '__main__':
	main()
