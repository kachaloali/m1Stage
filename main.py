# -*- coding: utf-8 -*-
import os
import sys
import glob
import commands
import operator
from optparse import OptionParser

parser = OptionParser(usage="%prog -f FILE")
parser.add_option("-f", "--files", dest = "files", help = "File name without fasta extension")
parser.add_option("-p", "--plma", dest = "plma", help = "PLMA's parameters", default= "-t 1 -c -M 7, -t 3 -c -M 7, -t 5 -c -M 7")
(args, options) = parser.parse_args()


# Paths to executables
pathprefix = os.path.dirname(sys.argv[0])
pathToClassificationExe = pathprefix +"/classification"
pathTo3DExe = pathprefix +"/3D"

#Files serving as a database
database1 = pathprefix +"/" + pathTo3DExe +"/database/pdbtosp.tab"
database2 = pathprefix +"/" + pathTo3DExe +"/database/uniprotData.tab"

#We retrieve the name of the fasta file
fileName = args.files

#The output file name
outputFile = fileName+'-res'

#The user is asked to pass the parameters of paloma if there are several, they must be separated by commas
#By default, if the user has not given any information, the parameters are: "-t 1 -c -M 7", "-t 3 -c -M 7", "-t 5 -c -M 7"
param = args.plma
if param != "" : param = param.split(',')
################################################################################################################################################
def main(fastaFile = fileName):
    
    if os.path.exists('./'+fileName+'_paloma'):
    	print commands.getoutput("rm -r "+fileName+'_paloma')
    print commands.getoutput("mkdir "+fileName+'_paloma')
    param2 = []
    
    #paloma
    print 'paloma ...'
    for i in range(len(param)):
        print param[i].strip()
	if os.path.exists(fastaFile + ".fasta"):
        	print commands.getoutput("paloma -i "+ fastaFile + ".fasta " + param[i])
        	ifile = glob.glob('./*.plma')[0]
	else: 
		print 'No such file or directory: '+ fastaFile + '.fasta'
		sys.exit(0)
        
       	liste, index1 = [], 0
        if "_" in ifile: liste = ifile.split("_")
        elem = [ elt for elt in liste if ".plma" in elt ]
        for elt in liste:
            if ".plma" not in elt: index1 += len(elt) + 1
            else: index2 = elt.find('.plma')
        index2 += index1
        param2.append(ifile[index1:index2])
            
        print commands.getoutput("mv ./*.plma "+fileName+'_paloma')
	print commands.getoutput("mv ./*.afc "+fileName+'_paloma')
	
    #classement
    classify(fastaFile, param, param2)
    
    #reclassement
    answer = raw_input("Do you want to reclassify [y/n]?\n")
    while answer in "oOyY":
        test = reclassify(fastaFile, param, param2)
        if test and answer not in "oOyY":
            answer = raw_input("Do you want to reclassify [y/n]?\n")
        else:
            print "There are no blocks to hide \n\n"
            answer = "n"
            
    #color3d.py
    answer = raw_input("Do you want to view PLMAs in 3D [y/n]?\n")
    if answer in "oOyY":
        print 'color 3d...'
        color3d(fastaFile, param, param2, database1, database2)
    else:
        print 'process terminated...'    
#=================================================================================================================================================
#function classify()
#=================================================================================================================================================
def classify(fastaFile, param, param2):
    #plma2dot
    print "dot ..."
    for i in range(len(param)):
        print param[i].strip()
        print commands.getoutput("plma2dot -i "+fastaFile+'_paloma/'+fastaFile+ "_"+ param2[i] +".plma")

    #colordot2.py
    for i in range(len(param)):
        print commands.getoutput("python ./"+ pathToClassificationExe +"/colordot2.py -d "+fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"_plma.dot -a "
                                                                +fastaFile+".attr -o " +fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"-col-out.dot")
    #plmadot2logic2.py
    print "logic ..."
    for i in range(len(param)):
	print param[i].strip()
        print commands.getoutput("python ./"+ pathToClassificationExe +"/plmadot2logic2.py -f "+fastaFile+".fasta -d "+fastaFile+'_paloma/'+fastaFile+"_"
                          +param2[i]+"-col-out.dot -a "+fastaFile+".attr > "+fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"-logic-out.lp")

    #blockclassconcept3.lp
    print 'concepts ...'
    for i in range(len(param)):
	print param[i].strip()
        print commands.getoutput("clingo "+ fastaFile +'_paloma/'+ fastaFile +"_"+ param2[i] +"-logic-out.lp ./"+ pathToClassificationExe
                                 +"/blockclassconcept3.lp --heu=Vsids --stats -n 0 > "+ fastaFile +"_paloma/trash ")
        print commands.getoutput("mv resblockclassconcept.lp "+ fastaFile +'_paloma/'+ fastaFile +"_"+ param2[i] +"-concepts-out.lp")

    #fullclassify5.lp
    print "classify ..."
    for i in range(len(param)):
	print param[i].strip()
        print commands.getoutput("clingo "+ fastaFile +'_paloma/'+ fastaFile +"_"+ param2[i] +"-concepts-out.lp ./"+ pathToClassificationExe
                                 +"/fullclassify5.lp --stats > "+ fastaFile +'_paloma/'+ fastaFile +"_"+ param2[i] +"-classif-out.lp")
    
    #ecrire2.py
    print "writing ..."
    command = "python ./"+ pathToClassificationExe +"/ecrire2.py -f "
    for i in range(len(param)):
        command += "./"+ fastaFile +'_paloma/'+ fastaFile +"_"+param2[i] + "-classif-out.lp,"
    command = command[:len(command)-1] + " -o " + outputFile +".csv -i "+ fastaFile +".fasta -a "+ fastaFile +".attr -p "+ pathToClassificationExe
    print commands.getoutput(command)
    if os.path.exists('./'+fileName+'_paloma/trash'):
        print commands.getoutput("rm ./"+ fastaFile +'_paloma/trash')
#=================================================================================================================================================
#function reclassify()
#=================================================================================================================================================
def reclassify(fastaFile, param, param2):
    #Fist, we write the blocks we want hidden in the classification file
    diUnusableBlocks = {}
    for i in range(len(param)):
        ifile = "./"+ fastaFile +'_paloma/'+ fastaFile +"_"+ param2[i] +"-classif-out.lp"
        xfile = open(ifile, 'r')
	lines = xfile.read().split('Answer:')[-1]
        if 'Optimization' in lines: lines = lines.split('Optimization')[0]
        elif 'Models' in lines: lines = lines.split('Models')[0]
        xfile.close()
        unusableBlocks = []
        blocks_new = list(filter(lambda line: 'blocknew(' in line, lines.split()))
        for line in blocks_new:
            numBlock = line.split(",")[1].strip(")")
            unusableBlocks.append(numBlock)
        diUnusableBlocks[i] = unusableBlocks
        ifile = "./"+fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"-concepts-out.lp"
        conceptfile = open(ifile, 'a')
        for block in unusableBlocks:
            conceptfile.write("hiddenblock(" + block +").\n")
        conceptfile.close()
        
    #fullclassify7.lp
    #Second, the classification is doing with some hidden blocks
    for i in range(len(param)):
        if diUnusableBlocks[i] != []:
            print commands.getoutput("clingo "+fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"-concepts-out.lp ./"+ pathToClassificationExe
                                     +"/fullclassify7.lp --configuration= handy --stats > "+fastaFile+'_paloma/'+fastaFile+"_"+param2[i]+"-classif-out.lp")

    #ecrire2.py         
    command = "python ./"+ pathToClassificationExe +"/ecrire2.py -f "
    for i in range(len(param)):
        if diUnusableBlocks[i] != []:
            command += "./"+fastaFile+'_paloma/'+ fastaFile +"_"+param2[i] + "-classif-out.lp,"
    command = command[:len(command)-1] + " -o " + outputFile +".csv -i "+ fastaFile +".fasta -a "+ fastaFile +".attr -p "+ pathToClassificationExe
    if "-classif-out.lp" in command:
        print commands.getoutput(command)
        return 1
    else:
        return 0
#=================================================================================================================================================
#function color3d()
#=================================================================================================================================================
def color3d(fileName, param, param2, database1, database2):
    if os.path.exists('./'+fileName+'_r3D'):
        print commands.getoutput("rm -r "+fileName+'_r3D')
    print commands.getoutput("mkdir "+fileName+'_r3D')

    print "searching for pdb structures ..."    
    print commands.getoutput("python ./"+ pathTo3DExe +"/sp2pdb.py -f "+ fileName +".fasta,"+ database1 +","+ database2 + " -p "+ pathTo3DExe)
        
    if os.path.exists('./obsolete'):
        print commands.getoutput("rm -r ./obsolete")
            
    for i in range(len(param)):
        print commands.getoutput("python ./"+ pathTo3DExe +"/color3d2.py -f "+ fileName +".fasta -d "+ fileName +'_paloma/'+ fileName +"_"+ param2[i]
                                                                                        +"-col.dot -p "+ fileName +"_pdb.attr --path "+ pathTo3DExe)
        print commands.getoutput("mv "+ fileName +"_paloma/*.pdf " + fileName +'_r3D')
        if os.path.exists("concepts"):
            print commands.getoutput("mv concepts " + fileName +'_r3D')
        
##################################################################################################################################################
	  
if __name__ == '__main__':
	main()




































