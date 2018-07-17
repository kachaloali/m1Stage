# -*- coding: utf-8 -*-
import os
import commands
import operator
from optparse import OptionParser

parser = OptionParser(usage="%prog -f FILE, FILE,... -o FILE -l FILE")
parser.add_option("-f", "--files", dest="files",help ="The classification files separated by commas")
parser.add_option("-o", "--out", dest="out",help ="The output file name")
parser.add_option("-i", "--fas", dest="fas",help ="The fasta file name")
parser.add_option("-a", "--attr", dest="attr",help ="The attibutes file PDB")
parser.add_option("-p", "--path", dest="path",help ="Path to programs")
(args, options) = parser.parse_args()

def main(files=args.files,output=args.out,fas=args.fas,attr=args.attr, pathToProg=args.path):
        #We retrieve the names of the classification files
        if ',' in files: files = files.split(',')
        else: files = files.split()
        
	diQueriesSeq, diNewFamily, param = {}, {}, []
	diAttrib, diNumSeqAndIdSeq = getIdSeqNumSeqAndColrs(fas,attr)
	fastaFileName = fas.replace('.fasta', '')
        if os.path.exists(fastaFileName +'_rClassif/'): print commands.getoutput('rm -r '+ fastaFileName +'_rClassif/')
	print commands.getoutput('mkdir '+ fastaFileName +'_rClassif/')
	####################################################################################################
	#We retrieve only the test sequences
	for idSeq, comment in diAttrib.items(): 
		if comment == 'black' :
                        diQueriesSeq[idSeq]=[]
                        for i in range(len(files)):
                                diQueriesSeq[idSeq].append([[], []])
        
        #For each file we replace each space with a line break and then retrieve the parameters of the file
        for ifile in files:
                print commands.getoutput("cat "+ ifile +" | tr \' \' \'\n\' > "+ ifile +'bis')
                print commands.getoutput("rm "+ ifile)
                print commands.getoutput("mv "+ ifile +'bis '+ ifile)
                #looking for the parameters
                liste, index1 = [], 0
                if "_" in ifile: liste = ifile.split("_")
                elem = [ elt for elt in liste if "-classif" in elt ]
                for elt in liste:
                        if "-classif" not in elt: index1 += len(elt) + 1
                        else: index2 = elt.find('-classif')
                index2 += index1
                param.append(ifile[index1:index2])
        ###################################################################################################              
	"""
        Here, if there are several classification files that are submitted, we run through each file and then recence
        the information provided. A sequence may be classified according to a classification file, and may not be 
	classified according to another file. It depends on the parameters used for the construction of these files.
        The parameters are those used since the alignment step (paloma)
        """
	diFile_concepts, counter = {}, 0
	for ifile in files:
                
                fileName, diBlocks, diTriFile, diClassement = ifile, {}, {}, {}
		xfile = open(ifile, 'r')
		if "/" in fileName:
                        chemin = fileName.split('/')[1]
                else:
                        chemin = os.getcwd()
		lines = xfile.read().split('Answer:')[-1]
		for iSeq in diQueriesSeq: diClassement[iSeq] = []
		#=========================================================================================
                if 'Optimization' in lines: lines = lines.split('Optimization')[0]; print 'Optimisation...'
                elif 'Models' in lines: lines = lines.split('Models')[0]; print 'Models...'
                #=========================================================================================
                bestclassified = list(filter(lambda line: 'bestclassified' in line.strip().split('(') and ',' in line.strip().split('(')[1], lines.split()))
                classified = list(filter(lambda line: 'classified' in line.strip().split('(') and ',' in line.strip().split('(')[1], lines.split()))
                bestambiguous = list(filter(lambda line: 'bestambiguous' in line.strip().split('(') and ',' in line.strip().split('(')[1], lines.split()))
                ambiguous = list(filter(lambda line: 'ambiguous' in line.strip().split('(') and ',' in line.strip().split('(')[1], lines.split()))
                unclassified = list(filter(lambda line: 'unclassified' in line.strip().split('('), lines.split()))
                new_family = list(filter(lambda line: 'support_new_family' in line, lines.split()))
                #=========================================================================================
                for line in bestclassified:
                        idSeq = (line.split(',')[0]).split('(')[1].strip('"')
                        if idSeq in diQueriesSeq:
                                diQueriesSeq[idSeq][counter][0].append(line.split(',')[1])
                                diQueriesSeq[idSeq][counter][1].append('best classified')
                                diClassement[idSeq].append(line.split(',')[1]) 
                                diTriFile[idSeq] = 6
               
                for line in classified:
                        idSeq = (line.split(',')[0]).split('(')[1].strip('"')
			if idSeq in diQueriesSeq:
                                diQueriesSeq[idSeq][counter][0].append(line.split(',')[1])
                                diQueriesSeq[idSeq][counter][1].append('classified')
                                diClassement[idSeq].append(line.split(',')[1])
                                diTriFile[idSeq] = 5

                for line in bestambiguous:
                        idSeq = (line.split(',')[0]).split('(')[1].strip('"')
                        if idSeq in diQueriesSeq:
                                diQueriesSeq[idSeq][counter][0].append(line.split(',')[1])
                                diQueriesSeq[idSeq][counter][1].append('best ambiguous')
                                diClassement[idSeq].append(line.split(',')[1])
                                diTriFile[idSeq] = 3


                for line in ambiguous:
                        idSeq = (line.split(',')[0]).split('(')[1].strip('"')
			if idSeq in diQueriesSeq:
                                diQueriesSeq[idSeq][counter][0].append(line.split(',')[1])
                                diQueriesSeq[idSeq][counter][1].append('ambiguous')
                                diClassement[idSeq].append(line.split(',')[1])
                                diTriFile[idSeq] = 2

                for line in unclassified:
                        idSeq = (line.split('("')[1]).strip('")')
                        if idSeq in diQueriesSeq:
                                diQueriesSeq[idSeq][counter][0].append('unclassified')
                                diQueriesSeq[idSeq][counter][1].append('')
                                diClassement[idSeq].append('unclassified')
                                diTriFile[idSeq] = 1
                ##################################################################################################                
                #Search for comcepts, associated blocks & associated sequences
                members_new = list(filter(lambda line: 'membernew(' in line, lines.split()))
                blocks_new = list(filter(lambda line: 'blocknew(' in line, lines.split()))
                test_quality = ['best classified', 'classified', 'best ambiguous', 'ambiguous', 'unclassified']
                diConcept = {}
                for line in new_family:
                        numConcept, iBlocks, iSeqs, infosConcept  = (line.split('(')[1]).split(',')[0], [], [], []

                        #The blocks members of the concept per file
                        blocks_of_concept = list(filter(lambda line: 'blocknew('+numConcept+',' in line,blocks_new))
                        for iline in blocks_of_concept:
                                numBlock = iline.split(',')[1].strip(')')
                                iBlocks.append(numBlock)
                        infosConcept.append(iBlocks)

                        #The sequences members of the concept per file
                        members_new_concept = list(filter(lambda line: ','+ numConcept +')' in line,  members_new))          
                        for iline in members_new_concept:
                                idSeq = iline.split('(')[1].split(',')[0].strip('"')
                                #If the sequence is among the queries sequences
                                if idSeq in diQueriesSeq:
                                        iSeqs.append(idSeq)
                                        diQueriesSeq[idSeq][counter][0].append('new('+ numConcept +')')
                                        if len(diQueriesSeq[idSeq][counter][1]) == 0:
                                                diClassement[idSeq].append('new('+ numConcept +')')
                                                diTriFile[idSeq] = 4
                                        
                        infosConcept.append(iSeqs)
                        diConcept[numConcept] = infosConcept
                diFile_concepts['File_'+str(counter+1)] = diConcept
                ##################################################################################################
                #Here we find the exceptions seauences ('except') if they exist.
                for idSeq in diQueriesSeq:
                        if len(diQueriesSeq[idSeq][counter][0]) == 0:
                                diQueriesSeq[idSeq][counter][0].append('except')
                                diClassement[idSeq].append('except')
                                diTriFile[idSeq] = 0
                                
                #Sorting the dictionary in descending order
                diTriFile = sorted(diTriFile.iteritems(), reverse=True, key=operator.itemgetter(1))
                if "/" in fileName:
                        outPutFile=open(fastaFileName+'_rClassif/'+fileName.split('/')[2].replace('classif-out.lp','res')+'.csv','w')
                else:
                        outPutFile=open(fastaFileName+'_rClassif/'+fileName.replace('classif-out.lp','res')+'.csv','w')
                outPutFile.write('File: '+fastaFileName+', param: '+ param[counter]+'\n\n\n')
                outPutFile.write('sequences   ,  subfamily       ,  quality  \n\n'.upper())
                
                #Writing results for each input classification file
                for i in range(len(diTriFile)):
                        idSeq = diTriFile[i][0]
                        outPutFile.write(idSeq+ ',')
                        for Class in list(set(diClassement[idSeq])) : outPutFile.write(Class + ' ')
                        outPutFile.write(','+ str(diTriFile[i][1]))
                        outPutFile.write('\n')

                
                xfileName = chemin+"/"+fastaFileName+"_"+param[counter]+"_plma.dot"
                diBlocks = getBlocks(xfileName)
                seqAndBlocks = getSeqAndInvolvedInBlocks(diNumSeqAndIdSeq,diBlocks)
                
                #Writing blocks
                outPutFile.write('\n\n  news families  \n\n\n'.upper())
                if diConcept != {}:
                        outPutFile.write("Concepts ,Members,Number of sequences,Number of blocks, interesting blocks\n")
                        for numConcept, conceptInfos in diConcept.iteritems():
                                if conceptInfos[1] !=[]:
                                        outPutFile.write(numConcept + ', ,'+ str(len(conceptInfos[1]))
                                                                 +','+ str(len(conceptInfos[0])) +'\n')
                                        for seq in list(set(conceptInfos[1])):
                                                suite_of_block = ''
                                                for numBlock in list(set(conceptInfos[0])):
                                                        if numBlock in seqAndBlocks[seq].keys():
                                                                suite_of_block += seqAndBlocks[seq][numBlock]+'  '
                                                outPutFile.write(","+ seq +',,,'+ suite_of_block+ "\n")
                        outPutFile.write('\n')
                        
                outPutFile.close()

                #Part Coloring PLMA by Families
                colorClassify(fas, attr, fileName, diQueriesSeq, diClassement, diConcept, param, counter, pathToProg)                
                
        	counter += 1
                xfile.close()
                
        """
	Writing step in the .csv file of the globals results, each sequence is written in the file with its status i.e
        Classified, ambiguous, unclassified etc. The subfamily field indicates the family (s) in which it was classified.
	"""
	outPutFile = open(fastaFileName+'_rClassif/'+output[:len(output)-4]+'Global'+output[len(output)-4:], 'w')
	outPutFile.write('File: '+fastaFileName+'\n\n\n')
	outPutFile.write('  sequences        , parameters    ,  subfamily       ,  quality  \n\n'.upper())

	for idSeq, infosSeq in diQueriesSeq.iteritems():
		outPutFile.write(idSeq)
		i = 0
		for liste in infosSeq:
                        outPutFile.write(',' + param[i] + ',')
			for Class in list(set(liste[0])) : outPutFile.write(Class + ' ')
			if len(liste[1]) > 0:
                                outPutFile.write(',' + liste[1][0] + '\n')
                        else:   outPutFile.write(', ' + '\n')
                        i +=1
 		outPutFile.write('\n')
	#For the new family
	outPutFile.write('\n\n  news families  \n\n\n'.upper())
	for File, Concept in diFile_concepts.iteritems():
                #=======================================================================================
                numFile = File[File.find('_')+1:]
                xfileName = chemin+"/"+fastaFileName+'_'+param[int(numFile)-1]+'_plma.dot'
                diBlocks = getBlocks(xfileName)
                seqAndBlocks = getSeqAndInvolvedInBlocks(diNumSeqAndIdSeq,diBlocks)
                #=======================================================================================
                if Concept != {}:
                        numFile = File[File.find('_')+1:]
                        outPutFile.write(File + ":  param  :  " + param[int(numFile) - 1]
                                         + ",Concepts ,Members,Number of sequences,Number of blocks, interesting blocks\n")
                        for numConcept, conceptInfos in Concept.iteritems() :
                                if conceptInfos[1] !=[]:
                                        outPutFile.write(','+ numConcept + ', ,'+ str(len(conceptInfos[1]))
                                                         +','+ str(len(conceptInfos[0])) +'\n')
                                        for seq in conceptInfos[1]:
                                                suite_of_block = ''
                                                for numBlock in list(set(conceptInfos[0])):
                                                        if numBlock in seqAndBlocks[seq].keys():
                                                                suite_of_block +=seqAndBlocks[seq][numBlock]+'  '
                                                outPutFile.write(", ,"+ seq +',,,'+ suite_of_block+ "\n")
                        outPutFile.write('\n')
	outPutFile.close()		
#########################################################################################################	
def getIdSeqNumSeqAndColrs(fas,attr):
    """
    This function returns two dictionaries where one of them, the keys are the id of the sequences & the values ​​are 
    the comments for each sequence. The other dictionary (diNumSeqAndIdSeq) its keys are the numbers of the sequences 
    in the PLMA file and the values ​​are the identifiers of the corresponding sequences.
    """
    with open(fas, 'r') as fFile:
        fastaFile=fFile.readlines()
        fFile.close()
    with open(attr, 'r') as aFile:
        attrFile=aFile.readlines()
        aFile.close()
    
    diQueriesSeq, diNumSeqAndIdSeq, numSeq = {}, {}, 0
    
    for fLine in fastaFile:
        if fLine[0] == '>':
            numSeq += 1

            if '|' in fLine:
                idSeq = fLine.split('|')[1].strip()
            else:
                idSeq = fLine[1:].strip()
            diQueriesSeq[idSeq] = ''
            diNumSeqAndIdSeq[str(numSeq)] = idSeq

            for aLine in attrFile:
                if 'range=' in aLine and 'comments=' in aLine:
                    borneInf = int(aLine.split('"')[1].split('-')[0])
                    borneSup = int(aLine.split('"')[1].split('-')[1])

                    if (borneInf <= numSeq and numSeq <= borneSup):
                        diQueriesSeq[idSeq] = aLine.split('"')[5]                  

    return diQueriesSeq, diNumSeqAndIdSeq
#################################################################################################
def getBlocks(dotFile):
	"""
	This function returns a dictionary of all the PLMA blocks contained in a dot file
        """
        with open(dotFile, 'r') as fd:
                dotfile = fd.readlines()
                
        subClustersDico = {}
        concatDotFile = reduce(lambda line1, line2: line1.strip()+line2.strip(), dotfile)
        subClusters = concatDotFile.split('subgraph cluster_')
  
        for subCluster in subClusters[3:]:
                subClusterTemp = subCluster.split('{')[1].split('"];')[:-1]
                tmp = subClusterTemp[0].strip().split(';')[2]
                subClusterTemp[0] = tmp
                subClustersDico[subCluster.split('{')[0]] = subClusterTemp
    
        lastSubCluster = subClusters[len(subClusters)-1:]
        lastSubClusterTemp = lastSubCluster[0].split('{')[1].split('}')[0].split('"];')[:-1]
        tmp = lastSubClusterTemp[0].strip().split(';')[2]
        lastSubClusterTemp[0] = tmp
        subClustersDico[lastSubCluster[0].split('{')[0]] = lastSubClusterTemp
     
        return subClustersDico
#################################################################################################
def getSeqAndInvolvedInBlocks(diNumSeq, diBlocks):
        diSeqBlocks = {}
        for numSeq, idSeq in diNumSeq.items():
                dico = {}
                for numblock, valueBlock in diBlocks.items():
                        for line in valueBlock:
                                if '"('+numSeq+', ' in line:
                                        dico[numblock] = line.split('label = "')[1]
                diSeqBlocks[idSeq] = dico
        return diSeqBlocks
##################################################################################################
def getNumSeqAndColrs(attribFile):
    """
    This function will make it possible to recover the sequence numbers and the color of their families
    """
    attributs = open(attribFile,'r')
    
    dico = {}
    for line in attributs.readlines():
        if 'range=' in line:
            ranger = line.split('"')[1]
            borneInf, borneSup = int(ranger.split('-')[0]), int(ranger.split('-')[1])
            color = line.split('"')[3]

            if borneInf > borneSup:

                error = "In the range section, the '-'  has to find "
                error += "between two numbers, and the first number "
                error += "has to be smaller than the second one!"
                printError(error)

            elif borneInf == borneSup:

                numSeq = borneInf
                dico[str(numSeq)] = color

            else:
                for numSeq in range(borneInf, borneSup+1):
                    dico[str(numSeq)] = color

    attributs.close()
    return dico

#################################################################################################
def colorClassify(fas, attr, fileName, diQueriesSeq, diClassement, diConcept, param, counter, pathToProg):
        fastaFileName = fastaFileName = fas.replace('.fasta', '')
        plma_seq1, plma_seq2  = getIdSeqNumSeqAndColrs(fas, attr)
        known_family = [family for family in list(set(plma_seq1.values())) if family != 'black']
        plma_seq3 = getNumSeqAndColrs(attr)
        colorNewFamily = "burlywood"
        colorAmbiguous = "olive"
        colorUnclassified = "black"

        diColor_of_family ={}
        for family in known_family:
                colors = []
                for numSeq in plma_seq3:
                        if plma_seq1[plma_seq2[numSeq]] == family.upper():
                                colors.append(plma_seq3[numSeq])
                diColor_of_family[family] = list(set(colors))

        
        colored_seq_by_family = {}
        for numSeq in plma_seq3:
                if plma_seq1[plma_seq2[numSeq]] != colorUnclassified:
                        colored_seq_by_family[numSeq] = []
                        colored_seq_by_family[numSeq].append(plma_seq3[numSeq])
             
        plma_seq2_temp = dict([[v,k] for v,k in plma_seq2.items()])
        
        #Inverting a dictionary               
        invert_dict = dict([[v,k] for k,v in plma_seq2.items()])
        plma_seq2 = invert_dict
        
        for idSeq in plma_seq1:
                if idSeq in diClassement:
                        numSeq = plma_seq2[idSeq]
                        colored_seq_by_family[numSeq] = []
                        for family, color_of_family in diColor_of_family.items():
                                if family.lower() in diClassement[idSeq]:
                                        colored_seq_by_family[numSeq].append(color_of_family[0])
                                
        colored_seq_by_family_tmp = dict([[cle,val] for cle,val in colored_seq_by_family.items()])
                                        
        #Give the color "colorNewFamily" for news families
        for idSeq in diClassement:
                for elem in diClassement[idSeq]:
                        if "new" in elem:
                                numSeq = plma_seq2[idSeq]
                                colored_seq_by_family[numSeq] = []
                                colored_seq_by_family[numSeq].append(colorNewFamily)
                       
        #Give the color "colorAmbiguous" for ambiguous
        for numSeq, list_color in colored_seq_by_family.items():
                if len(list_color) > 1:
                        colored_seq_by_family[numSeq] = []
                        colored_seq_by_family[numSeq].append(colorAmbiguous)

        
        #pools of family
        diFamily_by_colors = {}
        list_tmp = [ elem[0] for elem in colored_seq_by_family.values() if elem != [] ] 
        if colorNewFamily in set(list_tmp):
                diColor_of_family["new"] = [colorNewFamily]

        #Reverse of the dictionary of families and their colors
        invert_dict = dict([[v[0].lower(),k] for k,v in diColor_of_family.items()])
        diColor_family = invert_dict
        
        #A dictionary is created that contains the colors of the families and all the
        #sequences belonging to families
        for color_of_family in diColor_of_family.values():
                NumSeqs = []
                for numSeq, colorSeq in colored_seq_by_family.items():
                        if colorSeq != [] and colorSeq[0] == color_of_family[0]:
                                NumSeqs.append(numSeq)
                diFamily_by_colors[color_of_family[0]] = NumSeqs

        #Other unclassified sequences
        unclassified_seqs, list_tmp2 = [], []
        list_tmp1 = [ elem for elem in diFamily_by_colors.values()]
        for liste in list_tmp1:
                for elem in liste:
                        list_tmp2.append(elem)
                        
        list_tmp2 = list(set(list_tmp2))
        for numSeq in plma_seq3:
                if  numSeq not in list_tmp2:
                        unclassified_seqs.append(numSeq)
                        
        #Looking for ambiguous sequences
        ambiguous, reste_seqs, diClass = {}, {}, {}
        for numSeq, tColor in colored_seq_by_family.items():
                if numSeq in unclassified_seqs and tColor != []:
                        color = tColor[0]
                        ambiguous[numSeq] = color
                elif numSeq in unclassified_seqs:
                      reste_seqs[numSeq] = colorUnclassified
                      
        for numSeq in unclassified_seqs:
                color = colored_seq_by_family_tmp[numSeq] 
                if color != []: color = colored_seq_by_family_tmp[numSeq][0].lower()     
                else: color = ""
                if color != "":
                        if numSeq in colored_seq_by_family_tmp:
                                classes = diColor_family[color]
                                for color in colored_seq_by_family_tmp[numSeq][1:]:
                                        classes += ", " + diColor_family[color.lower()]
                        diClass[numSeq] = classes
       
        #==================================================================================================================
        #==================================================================================================================
        dotInFile = "./"+fastaFileName+"_paloma/"+fastaFileName+"_"+param[counter]+"_plma.dot"
        dotOutFile = "./"+fastaFileName+"_paloma/"+fastaFileName+"_"+param[counter]+"-col.dot"
        #==================================================================================================================
        #==================================================================================================================
        dic_blocks = {}
        lines = open(fileName, "r").readlines()
        #Looking for the characteristic blocks for each family
        for Class in diColor_of_family:
                blocks_support = list(filter(lambda line: 'characteristic_block' in line and Class.lower() in line, lines))
                blocks = []
                for line in blocks_support:
                        block = line.split(",")[2].split(")")[0]
                        blocks.append(block)
                dic_blocks[Class] = list(set(blocks))
                
        diChar_blocks = {}
        for Class, blocks in dic_blocks.items():
                for block in blocks:
                        diChar_blocks[block] = Class

        ####################################################################################################################
        #Creating of a dictionary that contains all the clusters of the plmadot
        dotFile = open(dotInFile, "r").readlines()
        subClustersDico, colorsSeq = {}, {}
        concatDotFile = reduce(lambda line1, line2: line1.strip()+line2.strip(), dotFile)
        subClusters = concatDotFile.split('subgraph cluster_')
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
        infoSeqs = lastSubCluster[0].split('{')[1].split('}')[1].split('];')[:-1]
        #===================================================================================================================
        #===================================================================================================================
        #The Input plmadot file
        inputFile = open(dotInFile, "r")
        #The output plmadot file
        outputFile = open(dotOutFile, "w")
        
        lines = inputFile.readlines()
        for index, elem in enumerate(lines):
                if "subgraph" in elem:
                        if elem.strip() == "subgraph cluster_1":
                                index1 = index
                        if elem.strip() == "subgraph cluster_2":
                                index2 = index
                        if elem.strip() == "subgraph cluster_3":
                                index3 = index
                                
        head = lines[:index1]                        
        cluster1 = lines[index1:index2]
        cluster2 = lines[index2:index3]
        
        #The sequences numbers and their labels
        diCluster1_tmp = {}
        for line in cluster1:
                if 'label' in line:
                        numSeq = line.split(",")[0].split('(')[1]
                        label = line.split(')"')[1]
                        diCluster1_tmp[numSeq] = label

        diCluster2_tmp = {}
        for line in cluster2:
                if 'label' in line:
                        numSeq = line.split(",")[0].split('(')[1]
                        diCluster2_tmp[numSeq] = line

        #===================================================================================================================
        #===================================================================================================================
        #The head of the dot is written
        for line in head:
             outputFile.write(line)
        #===================================================================================================================
        #===================================================================================================================
        #Part for cluster 1
        for line in cluster1:
                if "cluster" in line:
                        outputFile.write(line)
                        outputFile.write("{\n")
                elif "node" in line:
                        colorSeq = line.split('color =')[1].strip().split(',')[0]
                        line = line.replace(colorSeq.strip(), "black")
                        outputFile.write(line)
                elif "style" in line:
                        style_of_cluster = line.split("style =")[1].split(";")[0]
                        line = line.replace(style_of_cluster.strip(), "filled")
                        outputFile.write(line)
                        
        #Writing for the sub-families (cluster 1)
        i = 1
        allNewBlocks = []
        for color, NumSeqs in diFamily_by_colors.items():
                if color != colorNewFamily:
                        outputFile.write("subgraph cluster_" + str(i) +"p1 \n")
                        outputFile.write("{\n")
                        outputFile.write("label = \"Family: "+ diColor_family[color.lower()] +"\nNumber: "+ str(i) +"\";\n")
                        outputFile.write("node [shape = record, color =  black, fontcolor = black];\n")
                        for numSeq in NumSeqs:
                                if plma_seq2_temp[numSeq] in diQueriesSeq:
                                        line = diCluster1_tmp[numSeq].replace("\"];", " [**]\"];")
                                        outputFile.write('"('+numSeq+', 1, 0)"' + line)
                                else:   outputFile.write('"('+numSeq+', 1, 0)"' + diCluster1_tmp[numSeq])
                        outputFile.write('}\n')
                        i += 1
                #Case for pools of new families (if there are several)
                else:
                        i = 1
                        for concept, infosConcept in diConcept.iteritems():
                                outputFile.write("subgraph cluster_new" + str(i) +" \n")
                                outputFile.write("{\n")
                                outputFile.write("label = \"Family: "+ diColor_family[color.lower()] + "\nNumber: "+ str(i)
                                                 +"\";\n")
                                outputFile.write("node [shape = record, color =  black, fontcolor = black];\n")
                                for idSeq  in infosConcept[1]:
                                        numSeq = plma_seq2[idSeq]
                                        if idSeq in diQueriesSeq:
                                                line = diCluster1_tmp[numSeq].replace("\"];", " [**]\"];")
                                                outputFile.write('"('+numSeq+', 1, 0)"' + line)
                                        else:   outputFile.write('"('+numSeq+', 1, 0)"' + diCluster1_tmp[numSeq])
                                outputFile.write('}\n')
                                allNewBlocks += list(set(infosConcept[0]))
                                i += 1
        
        #We add the characteristic blocks of the new families
        for bloc in allNewBlocks:
                diChar_blocks[bloc] = "new"

        #The rest of the sequences (cluster 1)
        for line in cluster1:
                if 'label' in line: numSeq = line.split(",")[0].split('(')[1]
                if numSeq in unclassified_seqs:
                        color = colored_seq_by_family_tmp[numSeq]
                        if color != []:
                                color = colored_seq_by_family_tmp[numSeq][0].lower()
                        else: color = ""
                        if color != "":
                                if numSeq in colored_seq_by_family_tmp:
                                        classes = diColor_family[color]
                                        for color in colored_seq_by_family_tmp[numSeq][1:]:
                                                classes += ", " + diColor_family[color.lower()]
                                        line = line.replace(numSeq+ ':', "[" + classes.upper() +"] "+ numSeq+":")
                                        if plma_seq2_temp[numSeq] in diQueriesSeq:
                                                line = line.replace("\"];", " [**]\"];")
                                        outputFile.write(line)
                        else:
                                if plma_seq2_temp[numSeq] in diQueriesSeq:
                                        line = line.replace("\"];", " [**]\"];")
                                outputFile.write(line)            
        outputFile.write("}\n")
        #=================================================================================================================
        #=================================================================================================================
        #Part for cluster2
        for line in cluster2:
                if "cluster" in line:
                        outputFile.write(line)
                        outputFile.write("{\n")
                elif "node" in line:
                        colorSeq = line.split('color =')[1].strip().split(',')[0]
                        line = line.replace(colorSeq.strip(), "black")
                        outputFile.write(line)
                elif "style" in line:
                        style_of_cluster = line.split("style =")[1].split(";")[0]
                        line = line.replace(style_of_cluster.strip(), "filled")
                        outputFile.write(line)
                        outputFile.write("fontcolor = gray;\n")

        #Writing for the sub-families (cluster 2)
        i = 1
        for color, NumSeqs in diFamily_by_colors.items():
                if color != colorNewFamily:
                        outputFile.write("subgraph cluster_" + str(i) +"p2 \n")
                        outputFile.write("{\n")
                        outputFile.write("node [shape = record,style = filled, color =  "+color.lower()
                                                                                                +", fontcolor = black];\n")
                        outputFile.write("color = "+color.lower()+";\n")
                        for numSeq in NumSeqs:
                                outputFile.write(diCluster2_tmp[numSeq])
                        outputFile.write('}\n')
                        i += 1
                else:
                        i = 1
                        for concept, infosConcept in diConcept.iteritems():
                                outputFile.write("subgraph cluster_new" + str(i) +"\n")
                                outputFile.write("{\n")
                                outputFile.write("node [shape = record,style = filled, color =  "+color.lower()
                                                                                              +", fontcolor = black];\n")
                                outputFile.write("color = "+color.lower()+";\n")
                                for idSeq  in infosConcept[1]:
                                        numSeq = plma_seq2[idSeq]
                                        outputFile.write(diCluster2_tmp[numSeq])
                                outputFile.write('}\n')
                                i += 1
        
        #The rest of the sequences (cluster 2)
        for line in cluster2:
                if 'label' in line: numSeq = line.split(",")[0].split('(')[1]
                if numSeq in unclassified_seqs:  outputFile.write(line)
        outputFile.write("}\n")
        #=================================================================================================================
        #=================================================================================================================
        #Part for the rest of the clusters (PLMA blocks)
        for numCluster, cluster in subClustersDico.items():
                if numCluster in diChar_blocks:
                        outputFile.write("subgraph cluster_"+numCluster+"\n{\n")
                        outputFile.write("node [shape = record, style = filled, color = yellow, fontcolor = black];\n")
                        outputFile.write("color = "+diColor_of_family[diChar_blocks[numCluster]][0].lower()+";\n")
                        for line in cluster:
                                numSeq = line.split(",")[0].split("(")[1]
                                outputFile.write(line + "\"];\n")
                        outputFile.write("}\n") 
                elif numCluster not in ["1","2"]:
                        outputFile.write("subgraph cluster_"+numCluster+"\n{\n")
                        outputFile.write("node [shape = record, style = filled, color = yellow, fontcolor = black];\n")
                        outputFile.write("color = black;\n")
                        for line in cluster:
                                outputFile.write(line+"\"];\n")
                        outputFile.write("}\n")
             
        #Part for arrows
        for line in infoSeqs:
                if '->' in line:
                        numSeqTemp, numSeq = line.split('label = ')[1], ''
                        if ':' in line:
                                numSeq = numSeqTemp.split(':')[0].strip('"')
                        else:
                                numSeq = numSeqTemp.split(',')[0]
                        colorSeq = line.split(', color =')[1].strip().split(',')[0]
                        if numSeq in ambiguous:
                                line = line.replace("fontsize = 8","fontsize = 15")
                                line = line.replace("label = " + numSeq+ ',', "label = "+ numSeq +"("+ diClass[numSeq].upper()+")\"")
                                line = line.replace(colorSeq.strip(), ambiguous[numSeq].lower())
                                
                        elif numSeq in reste_seqs:
                                
                                color = plma_seq3[numSeq].lower()
                                if color != colorUnclassified:
                                        classe = diColor_family[color]
                                        line = line.replace("label = "+ numSeq+ ',', "label = \""+ numSeq+"("+ classe.upper() +")\"")
                                line = line.replace("fontsize = 8","fontsize = 15")
                                line = line.replace(colorSeq.strip(), "black") 
                        elif numSeq in colored_seq_by_family:
                                
                                if numSeq in colored_seq_by_family_tmp and colored_seq_by_family_tmp[numSeq] != []:
                                        
                                        color = plma_seq3[numSeq].lower()
                                        line = line.replace("fontsize = 8","fontsize = 15")
                                        if color != colorUnclassified:
                                                classe = diColor_family[color]
                                                line = line.replace("label = "+numSeq+ ',',"label = \""+ numSeq+" ("+ classe.upper() +")\"")
                                        else:
                                                line = line.replace("label = "+numSeq+ ',',"label = \""+ numSeq+" (?)\"")
                                                
                                elif colored_seq_by_family_tmp[numSeq] == []:
                                        color = colored_seq_by_family[numSeq][0]
                                        line = line.replace("fontsize = 8","fontsize = 15")
                                        classe = diColor_family[color]
                                        line = line.replace("label = "+numSeq+ ',',"label = \"" + numSeq+" (?)\"")
                                line = line.replace(colorSeq.strip(), colored_seq_by_family[numSeq][0].lower())
                                
                        outputFile.write(line+"];\n")
        outputFile.write("}\n")
        inputFile.close()
        outputFile.close()
        #================================================================================================================
        #================================================================================================================
        #Converting the product dot file to pdf format
        print commands.getoutput("python ./"+ pathToProg +"/plmadot2pdf.py -f ./"+fastaFileName+"_paloma/"+fastaFileName+"_"
                                                                                     + param[counter] +"-col.dot")
        print commands.getoutput("rm "+fastaFileName+"_paloma/"+fastaFileName+"_"+param[counter]+"-col.ps")
        print commands.getoutput("mv ./"+fastaFileName+"_paloma/"+fastaFileName+"_"+param[counter]+"-col.pdf ./"
                                                                                              +fastaFileName+"_rClassif")
#main
if __name__ == '__main__':
	main()
	
		
