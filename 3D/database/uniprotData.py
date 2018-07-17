# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:04:57 2017
@author: H. Kachalo Ali
"""

import ftplib
import subprocess
import ConfigParser


config = ConfigParser.ConfigParser()
config.readfp(open('conf.ini','r'))
uniprot_data = config.get('Download','UniProtKB')
UniprotIDTopdbId = config.get('Convert', 'UniprotIDTopdbId')
filename = "idmapping_selected.tab.gz"
#==================================================================================
#Remarque : le fichier est volumineux 4,3GB 
#==================================================================================
#On peut exécuter la commande suivante pour suivre l'évolution du téléchargement du
#fichier depuis la ligne de commande (le shell)
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
#==================================================================================
def downloadFile(filename):
    """
    - fonctionne pour la version 2.7 de python
    - se connecte au serveur ftp d'Uniprot et recupère le fichier compressé 
      au format tar.gz qui contient les informations que l'on souhaite a savoir: 
      les identifiant UniprotKB associés aux identifiants de la pdb
    """
    
    host = "ftp.uniprot.org" #adresse du serveur ftp
    connect = ftplib.FTP(host, 'anonymous', 'anonymous')
    connect.cwd('/pub/databases/uniprot/current_release/knowledgebase/idmapping/')

    with open(filename, 'wb') as infile:
        connect.retrbinary('RETR %s' % filename, infile.write)
        
    #décompression après avoir télécharger le fichier souhaité
    subprocess.check_output(['bash', '-c' , "gunzip -f " + filename])

def ExtractInfoFromFile(filename):
    """
    cette fonction permet d'extraire les informations qui nous concerne
    celles qui lient les identifiants uniprot au identifiants pdb
    """
    with open(filename.replace(".gz","") ,'r') as infile,\
    open(uniprot_data + ".tab", 'w') as output:
        uniprotID, PDBID = "", []
        for line in infile:
            lineList = line.split("\t")
            if (lineList[5] != "" and lineList[2] != ""):
                uniprotID = lineList[2]
                PDBIDList = lineList[5].split("; ")
                for pdbid in PDBIDList:
                    PDBID.append(pdbid)
                for pdbid in PDBID:
                    output.write(str(uniprotID) + "\t" + str(pdbid) + "\n" )
                uniprotID, PDBID = "", []

#appelle des fonctions
#==============================================================================
downloadFile(filename)             
ExtractInfoFromFile(filename)
#==============================================================================

