import commands
from optparse import OptionParser
parser = OptionParser(usage="%prog -f FILE")
parser.add_option("-f", "--dot", dest="dotf", metavar = "FILE" ,help ="fichier dot")
(args, options) = parser.parse_args()

def main(dotfile=args.dotf) :
    commands.getoutput("dot -Tps2 Gsize =\"100,100\" "+ dotfile +" -o "+ dotfile.replace('.dot', '') +".ps")
    commands.getoutput("ps2pdf "+ dotfile.replace('.dot', '') +".ps "+ dotfile.replace('.dot', '') +".pdf")

if __name__ == '__main__':
	main()
	
