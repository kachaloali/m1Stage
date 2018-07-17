# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 22:00:31 2017

@author: Kachalo ALi
"""
from myError import printError
from optparse import OptionParser
parser = OptionParser(usage="%prog -d FILE -c couleur1 -a FILE")
parser.add_option("-d", "--dot", dest="dotf", metavar = "FILE" ,help ="fichier dot")
parser.add_option("-a", "--attr", dest="attr", metavar = "FILE" ,help ="fichier attr")
parser.add_option("-f", "--fas", dest="fas", metavar = "FILE", help ="fichier fas")
(args, options) = parser.parse_args()

def main(dotF=args.dotf,attribfile = args.attr,opt=2) :
   
   with open(dotF, 'r') as dFile:
        dotFile=dFile.readlines()
        dFile.close()
    
   subClustersDico, colorsSeq = getSubclustersAndSeqColors(dotFile)
   seqsProtInfo = getProtNamesAndComments(subClustersDico, attribfile)
   write(subClustersDico, seqsProtInfo, colorsSeq)            
    
"""
permet de récuperer les blocks & les couleurs des séquences
"""
def getSubclustersAndSeqColors(dotFile):
    
    subClustersDico, colorsSeq = {}, {}
    concatDotFile = reduce(lambda line1, line2: line1.strip()+line2.strip(), dotFile)
    subClusters = concatDotFile.split('subgraph cluster_')
    
    #on recupère tous les blocks dans le dictionnaire (subClustersDico) dont les clés
    #sont les numeros des sous clusters du fichier plma.dot exp: subClustersDico['22']={..}
    for subCluster in subClusters[1:]:
        subClusterTemp = subCluster.split('{')[1].split('"];')[:-1]
        tmp = subClusterTemp[0].strip().split(';')[2]
        subClusterTemp[0] = tmp
        subClustersDico[subCluster.split('{')[0]] = subClusterTemp
    
    lastSubCluster = subClusters[len(subClusters)-1:]
    lastSubClusterTemp = lastSubCluster[0].split('{')[1].split('}')[0].split('"];')[:-1]
    tmp = lastSubClusterTemp[0].strip().split(';')[2]
    lastSubClusterTemp[0] = tmp
    subClustersDico[lastSubCluster[0].split('{')[0]] = lastSubClusterTemp
     
    #on recupère les numeros et les couleurs des séquences dans le dictionnaire (colorsSeq)
    #exp: colorsSeq[22]='green'
    infoSeqs = lastSubCluster[0].split('{')[1].split('}')[1].split('];')[:-1]
    for infoSeq in infoSeqs:
       
        numSeqTemp, numSeq = infoSeq.split('label = ')[1], ''
        colorSeq = infoSeq.split('color =')[1].strip().split(',')[0]
        
        if ':' in infoSeq:
            numSeq = numSeqTemp.split(':')[0].strip('"')
        else:
            numSeq = numSeqTemp.split(',')[0]
        colorsSeq[numSeq] = colorSeq 
        
    return subClustersDico, colorsSeq
   
"""
permet de récuperer les noms & classes(!si connues) des proteines 
"""
def getProtNamesAndComments(subClustersDico, attrFile):
    
    subClusterSeq = subClustersDico[str(1)]
    seqProtInfo = list()
    diComments = getNumSeqAndComments(attrFile)
    for val in subClusterSeq:
        
        if '|' in val:
            iSeq = val.split('label = "')[1].strip('"]').split('\|')
            seqInfos = list() 
            numSeq, chainSeq = iSeq[0].split(':')[0], '' 
            comments = diComments[str(numSeq)]
            
            if '_' in iSeq[1].strip(' '):
                idSeq = iSeq[1].strip(' ').split('_')[0]              
                chainSeq = iSeq[1].strip(' ').split('_')[1]
            elif ':' in iSeq[1].strip(' '):
                idSeq = iSeq[1].strip(' ').split(':')[0]              
                chainSeq = iSeq[1].strip(' ').split(':')[1]
            else:
                idSeq = iSeq[1].strip(' ')              
            
            seqInfos.append(idSeq.strip())#id de la sequence
            seqInfos.append(numSeq)#numero de la sequence dans le fichier plma.dot
            seqInfos.append(chainSeq)#la chaine de la sequence si'elle existe
            seqInfos.append(comments)#la classe de la sequence si'elle existe                
            seqProtInfo.append(seqInfos) 
            
        else:
            
            seqInfos = list() 
            iSeq = val.split('label = "')[1].strip('"]').split(':')
            numSeq = iSeq[0] 
            chainSeq, comments = '', diComments[str(numSeq)]
            idSeq = iSeq[1].strip(' ')
          
            seqInfos.append(idSeq.strip())#id de la sequence
            seqInfos.append(numSeq)#numero de la sequence dans le fichier plma.dot
            seqInfos.append(chainSeq)#la chaine de la sequence si'elle existe
            seqInfos.append(comments)#la classe de la sequence si'elle existe    
            seqProtInfo.append(seqInfos)               
    
    return seqProtInfo

'''
renvoie un dictionnaire dont les clés sont les num des séauences et les valeurs sont les comments  
'''
def getNumSeqAndComments(attribFile):
    attributs = open(attribFile,'r')
    dico = {}
    for line in attributs.readlines():
        if 'range=' in line:
            ranger = line.split('"')[1]
            borneInf, borneSup = int(ranger.split('-')[0]), int(ranger.split('-')[1])
            comments = line.split('"')[5].split('"')[0]

            if borneInf > borneSup:
                error = "Into the file " + attribFile + " in the range section, the '-'  has to find "
                error += "between two numbers, and the first number "
                error += "has to be smaller than the second one!"
                printError(error)

            elif borneInf == borneSup:
                numSeq = borneInf
                dico[str(numSeq)] = comments
            else:
                for numSeq in range(borneInf, borneSup + 1):
                    dico[str(numSeq)] = comments

    attributs.close()
    return dico                    

'''
permet d'ecrire les infos dans le fichier de sortie (séquences & blocs)
'''                    
def write(subClustDico, seqsProtInfo, colorsSeq):
   
    #écriture des séquences
    for infoSeq in seqsProtInfo:
        
        idSeq, numSeq, comments = infoSeq[0], infoSeq[1], infoSeq[3] 
        
        if comments == 'black':
            print 'sequence("' + str(idSeq) + '",black).'
        else:
            print 'sequence(' + str(numSeq) + ',' + comments.lower() + ').'
            
    #écriture des blocs
    for numCluster in range(3, len(subClustDico) + 1):
        for line in subClustDico[str(numCluster)]:
            
            numSeq = line.strip().split(',')[0].strip('"(')
            infoSeq = [line for line in seqsProtInfo if numSeq == line[1]]
            idSeq = infoSeq[0][0]
            comments = infoSeq[0][3]
            
            if comments == 'black':
                print 'block(' + str(numCluster) + ',"' + str(idSeq) + '").'
            else:
                print 'block(' + str(numCluster) + ',' + str(numSeq) + ').'                     
            



if __name__ == '__main__':
	main()
