#!/usr/bin/python
###############################################################################
# :: pfat.py ::
# Calvin Chen (cvchen@ucsd.edu)
# Joint Center for Structural Genomics
# @ San Diego Supercomputer Center, UC San Diego
#
# Written on GNU/Linux (RedHat 9.0)
#
#
# Align a FASTA peptide sequence with a PDB chain sequence. Format the output
# to reflect the true residue sequence positions of the PDB (sbjct) sequence,
# note gaps in the PDB sequence, note mismatches, and make the alignment
# "clickable" to display the sequence position of residues in the alignment.
# If necessary, "patch" the beginning and end of the alignment so that the
# entire FASTA sequence is displayed.
#
# There are two optional functions:
#	1) Create .pdb file where all residues in the original .pdb file are
#		renumbered according to the FASTA sequence it aligned with
#	2) Create .pdb file where the gaps in the PDB chain sequence is filled
#		with alanines.
###############################################################################

import os, re, sys
import tempfile
import commands
from optparse import OptionParser


blastPath = "/usr/local/bin/"
bl2seq_cmd = '%sbl2seq -p blastp -i %s -j %s -F F -X 50'

html = 'Content-type: text/html\n\n'
renum, gapfill, change = False, False, False

aas = { "ALA" : "A", "ARG" : "R", "ASN" : "N", "ASP" : "D", "CYS" : "C",
        "GLN" : "Q", "GLU" : "E", "GLY" : "G", "HIS" : "H", "ILE" : "I",
        "LEU" : "L", "LYS" : "K", "MET" : "M", "PHE" : "F", "PRO" : "P",
        "SER" : "S", "THR" : "T", "TRP" : "W", "TYR" : "Y", "VAL" : "V",
        "MSE" : "M"  }

chBlock = []

def alignpdbfas(fasFname, pdbFname, chainID = "A", outFname = "o.html", 
                rnFname = "", gfFname = ""):
   global blastPath, renum, gapfill, change
   outHTML = """<html>
<head>
<title>P-FAT Output</title>
<style = text/css">
<!--
A.aa:link{background-color:black; color:black; font-weight: none}
A.aa:hover{background-color:#666666; color:white; font-weight: bold}
BODY{
background-color:#CCCCCC;
font-family: Arial;
color: #666666;
}
PRE.down{
font-family: Arial;
font-weight: bold;
}
PRE.align{
background-color: #FFFFFF;
font-weight: none;
color: #000000;
}
-->
</style>
<script language="Javascript">
function show(query,subject){
   document.fields.que.value = query;
   document.fields.sub.value = subject;
}
</script>
</head>
<body>
<div align="center" style="font-family: Verdana; font-size: 46; color: #336699;"><b><u>P-FAT Output</u></b></div>
<u><b>Bold Underline</b></u>: Numerical gap in .pdb file between underlined letter and the one preceding it
<br><u><b><font color="#336699">Bold Italics Blue</font></b></u>: N and C terminal alignment extensions (not part of original alignment
<br><font color="red">Red</font>: Mismatch between sbjct and query ('-' does not count)
</form>
<hr>
<pre class="align"><!--align--></pre>
<form name="fields"><center>
<pre class="down">FASTA input (<i>Query</i>): <input name="que" type="text" size=10>		PDB input (<i>Sbjct</i>): <input name="sub" type="text" size=10></pre>
<!--renum--><!--gapfill-->
</center>
</form>
<hr>
Comments/questions: <a href="mailto:cvchen@ucsd.edu">Calvin Chen</a>, <a href="mailto:lukasz@sdsc.edu">Lukasz Jaroszewski</a>
<br><a href="http://www.jcsg.org" target="_blank">Joint Center for Structrual Genomics</a>
<br>@ <a href="http://www.sdsc.edu" target="_blank">San Diego Supercomputer Center</a> (<a href="http://www.ucsd.edu" target="_blank">UCSD</a>)
</body>
</html>"""


   # print '\nP-FAT: PDB-FASTA Alignment Tool'
   # print '\t--Calvin Chen\n\n'
   if not blastPath:
      print 'The Python and Blast path need to be set.'
      print 'Please run pfat_setup.py first.'

   # ask for and verify FASTA file
   if not os.path.isfile(fasFname):
      print '_Error_: %s is not a valid file.\n'%fasFname
   # ask for and verify PDB file
   if not os.path.isfile(pdbFname):
      print '_Error_: %s is not a valid file.\n'%pdbFname
   # ask for and verify chain ID
   if chainID.isalpha():
      chainID = chainID[0].capitalize()
   else:
      print '_Error_: %s is not a valid chain ID.\n'%chainID
   # ask for html output filename
   if not(outFname.endswith('.html') or outFname.endswith('.htm')):
       outFname = '%s.html'%outFname
   # option: re-number?
   if rnFname:
      renum = True
      change = True
      if not rnFname.endswith('.pdb'):
         rnFname = '%s.pdb'%rnFname
   # option: fill gaps?
   if gfFname:
      gapfill = True
      change = True
      if not gfFname.endswith('.pdb'):
         gfFname = '%s.pdb'%gfFname
     
   file = open(fasFname,'r')
   fasta_lines = file.read().split('\n')
   file.close()
   ffasta = ''
   if fasta_lines[0].startswith('>'):
      ffasta = ''.join(fasta_lines[1:])
   else:
      ffasta = ''.join(fasta_lines)
  
   file = open(pdbFname,'r')

   # getPDBseq() may get replaced by outside pdb-parsing script (ask Lukasz)
   cDict, seq, prev, seqNums = getPDBseq(file, chainID)
   if not seq:
      #print '%s does not contain a chain %s'%(pdbFname, chainID)
      sys.exit(0)
   else:
      # put PDB sequence into a file so we can BLAST it
      p_temp = doTemp('>%s Chain: %s\n%s'%(pdbFname, chainID, seq))

      # do alignment and format it
      align, pf = doAlign(
         bl2seq_cmd%(blastPath, fasFname, p_temp), seqNums, cDict, seq, ffasta)
      os.remove(p_temp)
      
      # there were options selected (re-numbering and/or gap filling)
      # process chBlock accordingly and write the new PDB files
      if change:
         theRest = file.read()
         if renum:
            newBlock = changePDB(pf)
            ofile = open(rnFname, 'w')
            ofile.write('%s%s%s'%(prev,newBlock,theRest))
            ofile.close()
         if gapfill:
            bridges = []
            # look for gaps to bridge
            for x in range(1,len(seqNums)):
               lo, hi = seqNums[x-1], seqNums[x]
               if hi - lo  > 1:
                  bridges.append((hi, bridgeGap(lo,hi,cDict)))
            if bridges:
               bridges.reverse()
               for x in range(0,len(chBlock)):
                  # patch the bridge into chBlock
                  if bridges and int(chBlock[x][22:26]) == bridges[-1][0]:
                     chBlock.insert(x, bridges.pop()[1])
            newBlock = ''.join(chBlock)
            ofile = open(gfFname,'w')
            ofile.write('%s%s%s'%(prev,newBlock,theRest))
            ofile.close()
      file.close()

      # make output html file
      sub = outHTML.replace('<!--align-->',align)
      file = open(outFname, 'w')
      file.write(sub)
      file.close()

###############################################################################
# When the pdb sequence is aligned with the fasta sequence via bl2seq, the
# residues that match may not have the same residue numbers. This function
# changes the block of text in the PDB file that describes the user specified
# chain and renumbers the amino acids so that each residue in that chain will
# have the same number as the amino acid that it matches with in the fasta
# sequence.
# -----------------------------------------------------------------------------
# Parameters:
#   nums -- (dict); maps pdb residue numbers to fasta residue numbers
#
# Return:
#   lines -- (str); the block of PDB text with changed amino acid numbers
###############################################################################
def changePDB(nums):
   lines = ''
   nums_keys = nums.keys()

   # each line in the block potentially needs to be modified
   for x in chBlock:
      num = x[22:26].strip()
      # we've got something to modify
      if num and int(num) in nums_keys:
         newNum = str(nums[int(num)]).rjust(4)
         lines = '%s%s%s%s'%(lines, x[:22], newNum, x[26:])
      # don't modify, just keep the line
      else:
         lines = '%s%s'%(lines,x)
   return lines

###############################################################################
# Use bl2seq on the fasta and pdb sequences, and format the output to display.
# Only the primary alignment is used; if there are secondary (etc) alignments,
# they are discarded. In the formatted alignment, the sequences are made to be
# 'clickable'; clicking on residues in the alignment will display the residues
# in three-letter representation and their sequence position. The residues of
# the PDB-derived sequence (Sbjct) are re-numbered to reflect their true
# numbering (as in the PDB file).
#
# Also, mismatches in the alignment are made red. Gaps in .pdb files exist
# when a residue in a chain has a position number that's more than one greater
# than the sequence position number of the residue preceding it. These gaps are
# denoted in the alignment by bold underlined sbjct sequence residues. The gap
# exists between the bold underlined residue and the one before it.
#
# bl2seq's behavior is to display the statistically significant
# portions of the two aligned sequences. Because this script is used to
# compare a (perhaps incomplete) protein in a PDB file against an assumably
# complete FASTA version, we want the displayed alignment to show the entire
# FASTA sequence, for comparison. The beginning and end of the FASTA sequence
# are patched into the alignment for this purpose, and they are either matched
# with unaligned beginnings/ends of the PDB sequence, or '-'s when they don't
# exist. No attempts are made to properly "align" these parts. These parts
# are distinguishable becaue the "Query" and "Sbjct" headings of the alignment
# lines are bold italics blue.
# -----------------------------------------------------------------------------
# Parameters:
#   cmd -- (str); the bl2seq command that aligns both sequences
#   aaNums -- (list); sequential list of amino acid numbers in .pdb file
#   cDict -- (dict); keyed by residue number, contains dictionaries of PDB file
#      lines further keyed by atom name
#   s_seq -- (str); the full pdb sequence
#   q_seq -- (str); the full fasta sequence
#
# Return:
#   (string); formatted output from bl2seq on pdb and fasta sequences
#   pf -- (dict); maps pdb residue numbers to fasta residue numbers
###############################################################################
def doAlign(cmd, aaNums, cDict, s_seq, q_seq):
   global chBlock
   query = re.compile(r'^Query:\s+(\d+)\s+([A-Z-]+)\s+(\d+)\s*$')
   sbjct = re.compile(r'^Sbjct:\s+(\d+)\s+([A-Z-]+)\s+(\d+)\s*$')
   identities = re.compile(r'^(\s+Identities = \d{1,3}/\d{1,3} \(\d{1,3}%\)).*$')
   atag = '<a class="aa" onClick="show(\'%s %s\',\'%s %s\')">%s</a>'
   aline = '%s <span style="cursor:default">%s</span> %s'
   block = ''
   sindex , qindex = 0,0
   pf = {}
   revAAs = {}
   s_gap = False
   s_2buf, q_2buf = [],[]
   startGap_s, startGap_q = 0,0
   bridges = []

   last_s, last_spos, last_qpos = '', -1000000, -1000000
   
   # invert the global aas dictionary into revAAs
   aas_keys = aas.keys()
   for key in aas_keys:
      revAAs[aas[key]] = key
   revAAs['X'] = 'ANY'

   out = commands.getoutput(cmd)

   # only keep parts of bl2seq output that we need
   sc = out.find('\n\n Score = ')
   if sc < 0:
    # print '\n\n__Error__: The selected .pdb chain does not have a significant'
    # print 'alignment with the FASTA sequence that you\'ve provided.'
    # print '\nbl2seq produced this output:\n\n:'
    # print out
     sys.exit(0)
   head = out[:out.find('\n\n Score = ')]
   outs = out[:out.find('\n\n\nLambda')].split('\n')

   # use only the primary alignment
   # alignments are separated by lines that contain "Score = "
   scount = 0
   cutTo = 0
   for x in range(0,len(outs)):
      theLine = outs[x]
      if theLine.find('Score = ') > 0:
         scount = scount + 1
         if scount > 1:
            outs = outs[:x]
            while(outs[-1].strip() == ''):
               outs = outs[:-1]
            break
         head = '%s\n\n%s'%(head,identities.match(outs[x+1]).group(1))
         cutTo = x+3
   outs = outs[cutTo:]

   very_last_qpos = int(outs[-3][outs[-3].rfind(' '):].strip())
   very_last_spos = int(outs[-1][outs[-1].rfind(' '):].strip())

   first_qpos = int(query.match(outs[0]).group(1))
   # fill in the missing beginnings of the seqeuences; do this by editing or
   # inserting to the front of outs[], it'll all get taken care of further down
   if first_qpos > 1:
      q_front = q_seq[:first_qpos-1]
      first_spos = int(sbjct.match(outs[2]).group(1))
      s_front = s_seq[:first_spos-1]
      s_temp = s_front[-len(q_front):]
      diffqs = len(q_front) - len(s_temp)
      s_front = '%s%s'%(diffqs*'-',s_temp)
      qposi, sposi = 1, first_spos - len(s_temp)
      toFront = []
      while(1):
         qline = q_front[:60]
         sline = s_front[:60]
         qstart = str(qposi).ljust(3)
         sstart = str(sposi).ljust(3)
         if len(qline) < 60:
            toFront.append('Query: %s %s %s'%(qstart, qline, first_qpos-1))
            toFront.append('')
            sEndNum = len(sline) - sline.count('-') + sposi
            toFront.append('Sbjct: %s %s %s'%(sstart, sline, sEndNum))
            toFront.append('')
            break
         qposi = qposi+60
         toFront.append('Query: %s %s %s'%(qstart, qline, qposi - 1))
         toFront.append('')
         sEndNum = len(sline) - sline.count('-') + sposi
         sposi = sEndNum + 1
         toFront.append('Sbjct: %s %s %s'%(sstart, sline, sEndNum))
         toFront.append('')
         q_front, s_front = q_front[60:], s_front[60:]
         if not q_front:
            break
      toFront.extend(outs)
      outs = toFront

   q_left = q_seq[very_last_qpos:]
   # fill in the missing ends of the sequences; do this by
   # editing/appending to outs[], it'll all get taken care of further down
   if q_left:
      qposTrack = very_last_qpos
      sposTrack = very_last_spos
      s_left = s_seq[sposTrack:]
      s_temp = s_left[:len(q_left)]
      diffqs = len(q_left) - len(s_temp)
      s_left = '%s%s'%(s_temp, diffqs*'-')

      while(1):
         new_last_q = q_left[:60]
         new_last_s = s_left[:60]
         last_q_start = str(qposTrack + 1).ljust(3)
         last_s_start = ''
         if new_last_s[0] == '-':
            last_s_start = str(sposTrack).ljust(3)
         else:
            last_s_start = str(sposTrack + 1).ljust(3)
         if len(new_last_q) < 60:
            outs.append('')
            outs.append('Query: %s %s %s'%(last_q_start, new_last_q, len(q_seq)))
            outs.append('')
            outs.append('Sbjct: %s %s %s'%(last_s_start, new_last_s, len(s_seq)))
            break
         qposTrack = qposTrack + 60
         sposTrack = sposTrack + 60 - new_last_s.count('-')
         outs.append('')
         outs.append('Query: %s %s %s'%(last_q_start, new_last_q, qposTrack))
         outs.append('')
         outs.append('Sbjct: %s %s %s'%(last_s_start, new_last_s, sposTrack))
         q_left, s_left = q_left[60:], s_left[60:]
         if not q_left:
            break

   # extract q/s line pairs, transform each line into 'clickable' html
   for x in range(0,len(outs)):
      match = query.match(outs[x])
      # found a query seq line, we know where the corresponding sbjct sequence
      # is; pull these two lines out and make them clickable
      if match:
         sseq, qseq = '', ''
         bgn = int(match.group(1))
         end = int(match.group(3))
         qindex = bgn
         q = match.group(2)                       # query
         match = sbjct.match(outs[x+2])
         s = match.group(2)                       # sbjct
         sbgn = int(match.group(1))
         sindex = sbgn-1
         
         # build clickable sequences for this pair of q/s lines
         for y in range(0,len(q)):
            # note starting aa positions for each row of each sequence
            spos = qpos = 0

            # get correct three leter aa format
            if q[y] == '-':
               qa = '---'
               qpos = qindex-1
            else:
               qa = revAAs[q[y]]
               qpos = qindex
               qindex = qindex + 1
            if s[y] == '-':
               sa = '---'
               if sindex-1 > 0 and sbgn < len(s_seq): 
                  spos = aaNums[sindex-1]
               else:
                  spos = aaNums[sindex]
               # note gap in PDB file
               if (spos - aaNums[sindex-2] > 1):
                  s_gap = True
            else:
               sa = revAAs[s[y]]
               spos = aaNums[sindex]
               # note gap in PDB file
               if (spos - aaNums[sindex-1] > 1):
                  s_gap = True
               sindex = sindex + 1
            if q[y] == '-':
               pf[spos] = 0
            elif s[y] is not '-':
               pf[spos] = qpos
            last_s, last_spos, last_qpos = s[y], spos, qpos
            # make this aa clickable and add onto growing sequence
            newqa = atag%(qa, qpos, sa, spos, q[y])
            newsa = atag%(qa, qpos, sa, spos, s[y])
            # we've got a mis-match; make it red
            if s[y] != '-' and q[y] != '-' and q[y] != s[y]:
               newsa = '<font color="red">%s</font>'%newsa
               newqa = '<font color="red">%s</font>'%newqa
            # denote gap with a bold underlined amino acid in alignment text
            if s_gap and s[y] != '-':
               newsa = '<b><u>%s</u></b>'%newsa
               s_gap = False
            qseq = '%s%s'%(qseq, newqa)
            sseq = '%s%s'%(sseq, newsa)

         ### add newly clickable sequences to block of clickable alignment text
         # this is part of the primary bl2seq alignment
         if bgn <= very_last_qpos and bgn >= first_qpos:
            block = '%s\n\n<b>Query</b>:%s\n<b>Sbjct</b>:%s'%(block,
               aline%(str(bgn).rjust(8),qseq,end),
               aline%(str(aaNums[sbgn-1]).rjust(8),sseq,spos))
         # this is not part of the primary bl2seq alignment
         else:
            font = '<font color="#336699"><i><u>%s</u></i></font>'
            block = '%s\n\n<b>%s</b>:%s\n<b>%s</b>:%s'%(block, font%'Query',
               aline%(str(bgn).rjust(8),qseq,end), font%'Sbjct',
               aline%(str(aaNums[sbgn-1]).rjust(8),sseq,spos))
      # look for 'Identities'
      else:
         match = identities.match(outs[x])
         if match:
            block = '%s\n\n\t%s'%(block,match.group(1))
   return '%s%s'%(head,block), pf

###############################################################################
# Given the higher and lower residue numbers in the sbjct sequence in a gap in
# the .pdb file, this function fills in the gap by introducing alanines with
# coordinates that in effect are a straight line connecting the two ends of the
# gap.
# -----------------------------------------------------------------------------
# Parameters:
#   slo -- (int); the lower sbjct residue number in the dash-gap
#   shi -- (int); the higher sbjct residue number in the dash-gap
#   cDict -- (dict); keyed by residue number, contains dictionaries of PDB file
#      lines further keyed by atom name
#
# Return:
#   bridge -- (str); chunk of lines that fill in particular gap
###############################################################################
def bridgeGap(slo, shi, cDict):
   tLine = 'ATOM  %5.5s  %sALA Z%4s    %8.3f%8.3f%8.3f  1.00 20.00           %.1s\n'
   ncacocb = ['N', 'CA', 'C', 'O', 'CB']
   atoms = [x for x in ncacocb if (x in cDict[slo] and x in cDict[shi])]
   dlo,dhi,diff = {},{},{}
   bridge = ''

   # compute dimensions and measurements of the gap to be bridged (filled)
   for a in atoms:
      loline = chBlock[cDict[slo][a]]
      hiline = chBlock[cDict[shi][a]]
      lox, loy, loz = loline[30:38].strip(), loline[38:46].strip(), loline[46:54].strip()
      hix, hiy, hiz = hiline[30:38].strip(), hiline[38:46].strip(), hiline[46:54].strip()
      lox, loy, loz = float(lox), float(loy), float(loz)
      hix, hiy, hiz = float(hix), float(hiy), float(hiz)
      dlo[a] = {'x': lox, 'y': loy, 'z': loz}
      dhi[a] = {'x': hix, 'y': hiy, 'z': hiz}
      diff[a] = {'x': hix - lox, 'y': hiy - loy, 'z': hiz - loz}

   gaplen = shi - slo

   # bridge the gap by computing atom coordinates
   for num in range(1,gaplen):
      for a in atoms:
         x = ((diff[a]['x']/gaplen)*num)+dlo[a]['x']
         y = ((diff[a]['y']/gaplen)*num)+dlo[a]['y']
         z = ((diff[a]['z']/gaplen)*num)+dlo[a]['z']
         bridge = '%s%s'%(bridge, tLine%(0,a.ljust(4),slo+num,x,y,z,a))
   return bridge

###############################################################################
# Make a temporary file and return it's file name for deletion later.
# -----------------------------------------------------------------------------
# Parameters:
#   write -- (str); contents of temprorary file
#
# Return:
#   temp_name -- (str); path and file name of newly created temp file
###############################################################################
def doTemp(write):
   tempfile.tempdir = os.getcwd()
   temp_name = tempfile.mktemp()
   temp = open(temp_name,'w')
   temp.write(write)
   temp.close()
   return temp_name

###############################################################################
# Get chain sequence from .pdb file. The lines in the .pdb file that come
# before the block of text for the chain we want is saved and returned, for use
# later in producing a modified .pdb file. 
# -----------------------------------------------------------------------------
# Parameters:  
#   pdbFile -- (str); file handle of the uploaded .pdb file
#   chainID -- (str); identify chain in .pdb file to align
#
# Return
#   cDict -- (dict); dict of PDB file lines of chainID, residue number maps to
#      dicts that map atom name to a specific PDB file line in chBlock
#   chainSeq -- (str); peptide sequence of chain from .pdb file
#   prev -- (str); text in pdb file preceding the chain we want to extract
#   aaNums -- (list); list of residue numbers in PDB chain
###############################################################################
def getPDBseq(pdbFile, chainID):
   flag = re.compile(r'^(HETATM|ATOM)\s*\d+.*$')
   end = re.compile(r'^TER\s+.*$')
   
   chainSeq = ''
   prev = ''
   cDict = {}
   global chBlock

   aaNum = -1000000  # this seems reasonably out of range
   aaNums = []
   aas_keys = aas.keys()

   # read through file and extract aa seqs of all chains
   while(1):
      line = pdbFile.readline()
      atom = line[12:16].strip()
      if not line:
         break
      # This is a file line that we are interested in looking at
      if(flag.match(line)):
         num = line[22:26].strip()
         # no amino acid number, probably passed the end of designated chain
         if chainSeq and not num:
            if change:
               # put this line back in file, it does not go with the block of
               # text that describes the chain we're after
               pdbFile.seek(-len(line),1)
            break
         newNum = int(num)
         # this is the chain we want
         if((line[21] == chainID) or (line[21] == 'A' and chainID == ' ')):
            # a valid amino acid 
            if(line[17:20] in aas_keys):
               if change:
                  if gapfill:
                     if cDict.has_key(newNum):
                        cDict[newNum][atom] = len(chBlock)
                     else:
                        cDict[newNum] = {atom: len(chBlock)}
                  chBlock.append(line)
               # new amino acid in our designated chain
               if (aaNum < newNum):
                  aaNum = newNum
                  aaNums.append(aaNum)
                  chainSeq = '%s%s'%(chainSeq,aas[line[17:20]])
            # invalid amino acid; probably passed the end of chain already
            else:
               if change:
                  # put this line back in file, it does not go with the block of
                  # text that describes the chain we're after
                  pdbFile.seek(-len(line),1)
               break
         # we're not past or even in the chain yet
         elif change:
            prev = '%s%s'%(prev,line)
      elif(end.match(line) and chainSeq and change):
         chBlock.append(line)
         break
      elif chainSeq:
         if change:
            # put this line back in file
            pdbFile.seek(-len(line),1)
         break
      elif change:
         prev = '%s%s'%(prev,line)
   return cDict, chainSeq, prev, aaNums



if __name__ == '__main__':
   parser = OptionParser(usage="%prog -i FILE -o FILE [options]")
   parser.add_option("-f", "--fas", dest="fasfile", metavar = "FILE", 
                     help = "fichier fas contenant la sequence")
   parser.add_option("-p", "--pdb", dest="pdbfile" , metavar = "FILE", 
                     help ="fichier pdb")
   parser.add_option("-c", "--chain", dest="chain",  default = 'A',
                     help = "chain ID")
   parser.add_option("-o", "--output", dest="outputfile", metavar = "FILE", 
                     default = 'o.html' , help = "output html file")
   parser.add_option("-O", "--outputpdb", dest = "outputpdbfile",
                     metavar = "FILE", help = "output pdb file")
   parser.add_option("-g", "--gapfill", dest = "gapopt", default = '', 
                     help ="option gapfill")
   (args, options) = parser.parse_args()
   alignpdbfas(args.fasfile, args.pdbfile, args.chain, args.outputfile, 
               args.outputpdbfile, args.gapopt)

