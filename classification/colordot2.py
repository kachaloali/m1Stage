import argparse

parser = argparse.ArgumentParser(prog = "colorDot-Ali.py")
parser.add_argument("-d", "--dot", dest = "dotf", metavar = "FILE", help = "dot file", required = True)
parser.add_argument("-a", "--attr", dest = "attr", metavar = "FILE", help = "attribute file", required = True)
parser.add_argument("-o", "--output", dest = "outf", metavar = "FILE", help = "output file", default = "colorDot-out.dot")
parser.add_argument("--version", action = "version", version = "%(prog)s 1.0")
args = parser.parse_args()

def main(dotInFile = args.dotf, attrFile = args.attr, dotOutFile = args.outf):

    inputFile = open(dotInFile, "r")
    outputFile = open(dotOutFile, "w")
  
    diColors = getNumSeqAndColrs(attrFile)

    for line in inputFile.readlines():
        if '->' in line:

            numSeqTemp, numSeq = line.split('label = ')[1], ''

            if ':' in line:
                numSeq = numSeqTemp.split(':')[0].strip('"')
            else:
                numSeq = numSeqTemp.split(',')[0]

            colorSeq = line.split('color =')[1].strip().split(',')[0]

            line = line.replace(colorSeq.strip(), diColors[numSeq])

        outputFile.write(line)

    inputFile.close()
    outputFile.close()
    

def getNumSeqAndColrs(attribFile):

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




if __name__ == '__main__':
	main()
