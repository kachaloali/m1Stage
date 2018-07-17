# -*- coding: utf-8 -*-
"""
Created on Wed May 24 08:40:31 2017

@author: H. Kachalo Ali
This program needs pdbtosp.txt.html in current directory
(to get/update it you may need to run the command: wget -c http://www.uniprot.org/docs/pdbtosp.txt.html to download pdbtosp.txt.html)

After computed this program, the output file is: 'pdbtosp.tab'

Number of PDB entries referenced in Swiss-Prot: 95913
Number of Swiss-Prot entries with one or more pointers to PDB: 24786
"""

with open('pdbtosp.txt.html','r') as fd:
    
    infile = fd.readlines()
    outputfile = open('pdbtosp.tab', 'w')
    outputfile.write("\n\nNumber of PDB entries referenced in Swiss-Prot: 95913\n")
    outputfile.write("Number of Swiss-Prot entries with one or more pointers to PDB: 24786\n\n\n")
    outputfile.write('pdbId\tMethod\tResolution\tsp entry name(s)\tspId\n')
    i = 0
    while i < len(infile):
        line = infile[i]
        if "pdbe-srv/view/entry/" in line:
            iline = line.split('">')
            if len(iline) == 3:
                
                iline1, iline2 = iline[1], iline[2]
                pdbcode = iline1.split('<')[0]
                method = iline1.split()[1]
                resolution = iline1.split()[2]
                if  resolution == '-':  sProtEntry = iline1.split()[3]
                else:  sProtEntry = iline1.split()[4]
                idProtein = iline2.split('<')[0]
                outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                 +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
            else:
                
                iline = line.split(',')
                iline0 = iline[0].split('">')
                iline1, iline2 = iline0[1], iline0[2]
                pdbcode = iline1.split('<')[0]
                method = iline1.split()[1]
                resolution = iline1.split()[2]
                if  resolution == '-':  sProtEntry = iline1.split()[3]
                else:  sProtEntry = iline1.split()[4]
                idProtein = iline2.split('<')[0]
                outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                 +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                iline_rest = iline[1:]
                for prot_info in iline_rest:
                    
                    if len(prot_info) > 1:
                        idProtein = prot_info.split('">')[1].split('<')[0]
                        sProtEntry = prot_info.split()[0].strip()
                        outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                     +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                    else:
                        prot_info = infile[i+1]
                        idProtein = prot_info.split('">')[1].split('<')[0]
                        sProtEntry = prot_info.split()[0].strip()
                        outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                     +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                        i +=1
                    
        else:
            
            if('http://www.uniprot.org/uniprot/') in line:
                iline = line.split('">')
                if len(iline) == 2:
                    
                    iline1, iline2 = iline[0], iline[1]
                    pdbcode = iline1.split()[0]
                    method = iline1.split()[1]
                    resolution = iline1.split()[2]
                    if  resolution == '-':  sProtEntry = iline1.split()[3]
                    elif 'X-ray' in iline1:  sProtEntry = iline1.split()[4]
                    idProtein = iline2.split('<')[0]
                    outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                     +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                else:
                    
                    iline = line.split(',')
                    iline0 = iline[0].split('">')
                    iline1, iline2 = iline0[0], iline0[1]
                    pdbcode = iline1.split()[0]
                    method = iline1.split()[1]
                    resolution = iline1.split()[2]
                    if  resolution == '-':  sProtEntry = iline1.split()[3]
                    elif 'X-ray' in iline1:  sProtEntry = iline1.split()[4]
                    idProtein = iline2.split('<')[0]
                    outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                     +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                    iline_rest = iline[1:]
                    for prot_info in iline_rest:
                        
                        if len(prot_info) > 1:
                            idProtein = prot_info.split('">')[1].split('<')[0]
                            sProtEntry = prot_info.split()[0].strip()
                            outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                             +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                        else:
                            prot_info = infile[i+1]
                            idProtein = prot_info.split('">')[1].split('<')[0]
                            sProtEntry = prot_info.split()[0].strip()
                            outputfile.write(pdbcode +'\t'+ method +'\t'+ resolution
                                             +'A\t'+ sProtEntry +'\t'+idProtein+'\n')
                            i +=1
        i +=1
    outputfile.close()
    print 'the end.'
#############################################################################################################
