from Bio.PDB.PDBParser import PDBParser
import numpy
#import mixture, random
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
#import fastcluster
import scipy.signal as signal
import math 

#import sklearn.cluster

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

import copy
#import swalign

from scipy.fftpack import *

import scipy
import argparse
import os
#################################################

CONTACT_TS = 7.5
INTERACTION_TS = 11
padSize = 46

homeDir = '/home/cgaliez/'
rootDir = homeDir + 'Tests/'
PDBDirectory = homeDir + 'Downloads/PDB/'
HBondsDirectory = rootDir + 'Results/HBonds/'
MIDirectory = rootDir + 'Results/MutualInformation/'
InterpolatedMotifDirectory = rootDir + 'Results/Motifs/Interpolated/'
PaddedMotifDirectory = rootDir + 'Results/Motifs/Padded/'
MotifDirectory = rootDir + 'Results/Motifs/'
SequenceMotifDirectory = rootDir + 'Results/SeqMotifs/'
centroidsSeqDirectory = rootDir + 'Results/Motifs/CentroidsSeqModels/'
MotifVizDirectory = rootDir + 'Results/MotifsViz/'
MotifVizFFTDirectory = rootDir + 'Results/MotifsVizFFT/'
MotifsBySizeDirectory = rootDir + 'Results/Motifs/MotifsBySize/'
SeqFASTAMotifsDirectory = rootDir + 'Results/SeqFASTAMotifs/'
SeqFASTAMotifsFFTDirectory = rootDir + 'Results/SeqFASTAMotifsFFT/'
SeqMotifsDirectory = rootDir + 'Results/SeqMotifs/'
MotifModelsDir = rootDir + 'Results/MotifsModels/'
someMotifFFTDirectory = rootDir + 'Results/Motifs/MotifsFFT/'

#mode  = 'interpolated'
mode  = 'padded'

if mode == 'padded':
	someMotifDirectory = PaddedMotifDirectory
	imageMotifDirectory = MotifDirectory + 'PaddedMotifsImages/'
elif mode == 'interpolated':
	someMotifDirectory = InterpolatedMotifDirectory
	imageMotifDirectory = MotifDirectory + 'MotifsImages/'

#################################################


# List of amino acids with three letter code.
AAIndex = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

AANick = { 'ASP' : 'D' , 'GLU' : 'E' , 'GLN' : 'Q' , 'ASN' : 'N' , 'SER' : 'S' ,
       'THR' : 'T' , 'CYS' : 'C' , 'HIS' : 'H' , 'ARG' : 'R' , 'LYS' : 'K' ,
       'MET' : 'M' , 'ALA' : 'A' , 'ILE' : 'I' , 'LEU' : 'L' , 'VAL' : 'V' ,
       'GLY' : 'G' , 'PRO' : 'P' , 'TRP' : 'W' , 'PHE' : 'F' , 'TYR' : 'Y' ,
       'SCY' : 'U' , 'ASX' : 'B' , 'GLX' : 'Z' , 'XXX' : 'X'}

AANickIndex = [AANick[x] for x in AAIndex]

# The last AA is a virtual AA placed at the centre of mass of all AA (may
# should be ponderated by the pb of occ of each type of AA
AAVect = numpy.array([[-0.2178,0.7889,-1.2688,-1.1404,1.3707,-0.6713,-0.5473,1.2621,-0.9216,1.1726,0.9232,-1.0798,-0.5130,-0.8079,-0.7383,-0.5972,-0.1661,0.9925,1.2249,0.9344],
[0.6688,0.9930,-0.0232,-0.1138,-0.7666,0.1578,-1.2124,0.8360,-0.0573,0.4977,0.4925,-0.0821,0.4129,-0.3783,-0.3316,0.4368,0.6951,0.9202,-1.6356,-1.5099],
[-0.4686,-1.3339,-0.3380,0.4484,-0.1158,-1.2349,0.2896,0.5295,0.7535,0.6694,1.0309,-0.3552,-0.0784,0.8135,0.7914,-0.5908,-0.1994,0.3519,-0.7134,-0.2497],
[-0.3641,0.6450,0.4331,-0.0896,0.3056,-0.6443,0.7553,0.2185,-0.1954,-0.0059,-0.1208,0.8978,-1.6099,-0.0648,-0.1989,0.1649,0.4102,0.0293,-1.0520,0.4859],
[-0.4557,-0.7067,1.3794,0.7423,0.3287,-0.0818,-0.0700,0.4381,-0.6943,0.2635,-0.0941,0.1004,0.6359,-0.3766,-1.1079,-0.3439,-0.1771,0.2605,-0.1610,0.1204],
[0.1424,0.4257,-0.4105,-0.3337,0.2818,-0.9743,0.7453,-0.0863,-0.0324,-0.2390,-0.5159,-0.2698,1.4979,-0.2103,-0.0125,0.1851,0.1549,-0.1600,-0.8044,0.6159],
[-0.0344,0.8548,0.2437,0.9006,-0.4543,-1.2470,-0.4752,-0.0480,0.0125,-0.4385,-0.3620,-0.3046,-0.1311,0.5916,-0.1294,-0.0101,0.1181,0.2485,0.7996,-0.1348],
[0.6992,-1.0898,-0.5216,0.1185,-0.1822,-0.2027,0.1805,-0.0432,-0.4504,-0.2557,-0.0472,-0.1553,-0.2878,0.1631,-0.6192,0.9250,1.0355,0.3534,0.2027,0.1770],
[-0.0360,-0.4101,0.3228,0.0327,0.6636,-0.2061,-1.0415,0.1057,0.7799,-0.1480,-0.6666,-0.1680,-0.1828,-0.7158,0.6287,0.2494,0.3906,0.1185,-0.1341,0.4173],
[0.7730,0.1321,-0.0561,0.6490,0.2578,0.4701,0.2891,0.0633,0.0976,0.0376,-0.4613,-1.0377,-0.4624,0.2323,-0.1341,-0.1000,-0.6941,0.2388,-0.5510,0.2559],
[-0.1088,-0.1509,-0.1395,-0.0102,-0.7825,0.1482,0.6186,0.7195,-0.0818,-0.2172,-0.7016,0.1668,-0.0183,-0.6130,0.5061,-0.3895,-0.0287,0.8283,0.3024,-0.0477],
[0.1385,-0.1176,-0.3490,-0.2885,-0.0052,0.0993,-0.4225,0.2601,0.4897,-0.2382,-0.3316,0.6900,0.1201,0.6577,-0.6273,-0.3456,-0.3024,0.3303,-0.1165,0.3587],
[-0.7961,0.1683,-0.2264,0.1878,0.0572,0.5346,0.0710,0.0661,0.1607,-0.1517,-0.1339,-0.4827,0.0482,0.3273,-0.1865,-0.4055,0.8832,-0.0468,-0.1476,0.0730],
[0.1699,0.0616,0.1294,0.0672,-0.1180,0.0233,0.0673,-0.0272,0.3286,-0.9159,0.8264,-0.1150,-0.0159,-0.4367,-0.1772,-0.2283,-0.0077,0.1407,-0.0501,0.2777],
[0.5807,-0.0498,0.5113,-0.4940,-0.1731,-0.1016,-0.0419,-0.2547,-0.1538,0.1738,-0.0552,-0.1522,-0.0444,0.2130,0.1611,-0.7320,0.4404,0.0038,0.0092,0.1595],
[-0.0906,-0.0456,-0.0349,0.0234,0.2955,0.1238,-0.3036,0.0963,-0.8188,-0.5338,0.0917,0.0960,0.0818,0.4017,0.5841,0.0082,-0.0984,0.2399,-0.2183,0.1016],
[0.2763,0.0143,-0.4482,0.5344,0.4914,-0.0540,0.0636,-0.3953,-0.0690,0.0462,-0.0417,0.4366,0.0283,-0.3075,0.0066,-0.5418,0.1998,0.1172,0.0124,-0.3696],
[-0.0244,-0.0480,0.2719,-0.4487,0.5703,-0.0900,0.3356,0.0654,0.2181,-0.3573,-0.1338,-0.1735,-0.0315,0.1954,-0.1121,0.0877,-0.0299,0.2765,0.1618,-0.7333],
[0.2786,-0.0303,-0.0596,0.0890,0.0911,-0.0392,0.0031,0.7120,-0.0443,-0.1745,-0.0447,0.0427,-0.0243,0.0081,0.0160,-0.1455,0.1027,-0.7048,0.0687,-0.1447],
[0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000],
]).T
#[0.11309,0.0101,-0.05842,0.08744,0.2116,-0.39901,-0.06954,0.45179,-0.06785,-0.08149,-0.03457,-0.19456,-0.05747,-0.03072,-0.1681,-0.23731,0.27267,0.45387,-0.28023,0.07876]


# Some caracteritic poles in the AA space, here, (currently) are hydrophobic and hydrophilic poles
AAPoles = [[0.1792111111111111, -0.012411111111111126, -0.046388888888888875, -0.012699999999999982, 0.098588888888888906, -0.034322222222222219, -0.044755555555555548, 0.45035555555555551, -0.21488888888888888, 0.10741111111111112, -0.0074444444444444384, -0.24650000000000002, -0.032055555555555552, -0.17521111111111112, -0.28012222222222227, -0.014066666666666647, 0.081933333333333344, 0.35171111111111114, -0.043011111111111111, -0.10532222222222221]
,[-0.029459999999999997, 0.089120000000000005, 0.10047999999999999, 0.067719999999999975, 0.017639999999999996, -0.18416000000000002, -0.38202000000000003, 0.025840000000000009, 0.25780000000000003, -0.16052, -0.30985999999999997, -0.083500000000000019, -0.037999999999999992, 0.21475999999999998, -0.030679999999999978, -0.24876000000000001, 0.30598000000000003, 0.13086000000000003, 0.082119999999999999, 0.17473999999999998]
]

#############################
# Secondary structure class #
#############################

class SecStruct:
	def __init__(self,t,s,e,d=0):
		self.ssType = t
		self.start = s
		self.end = e
		self.direction = d

	def __repr__(self):
		return repr((self.ssType,self.start,self.end,self.direction))

	def __str__(self):
		return {'H' : 'Helix','S' : 'Strand','C': 'Coil'}[self.ssType] + ' from ' + str(self.start) + ' to ' + str(self.end) + {0 : '' , 1 : ' in parallel sense', -1 : ' in anti-parallel sense'}[self.direction]





DEFAULT_THRESHOLD = 7.5
#############################
# Conatct map class         #
#############################
class contactMap:
	def __init__(self,arg,thsld = DEFAULT_THRESHOLD):

		self.seqCutOff = 4
		if type(arg) is str: # arg is cm file name
			self.parseCM(arg)
		elif arg.shape[0] == arg.shape[1]: # distance matrix
			self.createCMFromDM(arg,thsld)
		else: # array of coords
			self.createCM(arg,thsld)

	def createCM(self,coords,thsld):

		# Create contact map
		self.cm = list()
		self.ts = thsld
		for i in range(len(coords)):
			for j in range(i,len(coords)):
				if abs(i-j) >= self.seqCutOff and ((coords[i][0] - coords[j][0])**2 + (coords[i][1] - coords[j][1])**2 + (coords[i][2] - coords[j][2])**2)**.5 < self.ts:
					self.cm += [(i,j)]
		self.cm = numpy.array(self.cm)


	def createCMFromDM(self,dm,thsld):
		# Create contact map
		self.cm = list()
		self.ts = thsld
		N = dm.shape[0]
		for i in range(N):
			for j in range(i,N):
				if abs(i-j) > self.seqCutOff and dm[i,j] < self.ts:
					self.cm += [(i,j)]
		self.cm = numpy.array(self.cm)

	def parseCM(self,filename):
		f = open(filename, 'r')

		self.cm = list()

		for line in f.readlines():
			if line[0] != '#':
				self.cm += [(int(line.split()[0])-1,int(line.split()[1])-1)]
		f.close()
		self.cm = numpy.array(self.cm)


	def __repr__(self):
		return self.cm

	def __str__(self):
		out = ''
		for (i,j) in self.cm:
			out += str(i)+ " --> " + str(j) + '\n'
		return out



#############################
# PDB parser                #
#############################

class myPDBParser:
	def __init__(self,filename,chainsToParse='all'):
		self.hbondAngleTS = 90
		self.hbondTS = 3.5
		self.dmCreated = False
		self.HBonds = list()
		self.HBondsParsed = False

		if filename != '':
			self.pdbCode = filename.split('/')[-1].split('.')[0] # if the filename has the format '/path/to/pdb/code.pdb'
			self.filename = filename
			self.parseCoords(filename,chainsToParse)
			#self.sort()


	def parseSS(self,filename):
		f = open(filename, 'r')
		self.ss = []

		for line in f.readlines():
			if line.split()[0] == 'HELIX':
				start,end = int(line.split()[5]),int(line.split()[8])
				self.ss += [SecStruct('H',start,end)]
			if line.split()[0] == 'SHEET':
				start,end,direction = int(line.split()[6]),int(line.split()[9]),int(line.split()[10])
				self.ss += [SecStruct('S',start,end,direction)]

	# Rectangle is given in residue numbering mode (i.e. index written in the PDB, *not* the index in our internal res list)
	def getPDBExtract(self,rectangle):
		f = open(self.filename, 'r')
		(xm,chainxm),(xM,chainxM),(ym,chainym),(yM,chainyM) = rectangle
		extract = list()
		lineCount = 0
		for line in f.readlines():
			fields = line.split()
			if len(fields) > 0 and fields[0] == 'ATOM':
				# if the line describe a carbon alpha and in the rectangle
				if fields[2] == 'CA':
					if ((fields[4]==chainxm or fields[4]==chainxM) and int(fields[5]) >= xm and int(fields[5]) <= xM) or ((fields[4] == chainym or fields[4]==chainyM) and int(fields[5]) >= ym and int(fields[5]) <= yM):
						extract+= [line] # then store it
					lineCount+=1
		f.close()
		return ''.join(extract)

	# Rectangle is given in residue numbering mode (i.e. index written in the PDB, *not* the index in our internal res list)
	def writePDBExtract(self,rectangle,outFile):
		f = open(outFile,'w')
		f.write(self.getPDBExtract(rectangle))
		f.close()

	def parseCoords(self,filename, chainsToParse='all'):
		p=PDBParser(PERMISSIVE=1,QUIET=True)
		s=p.get_structure("", filename)

		self.coords = list()
		self.residues = list()
		self.sequence = list()
		self.resID = list()
		self.elements = dict()
		self.chainEntries = dict()
		self.chainLen = dict()
		self.elements["N"] = list()
		self.elements["F"] = list()
		self.elements["S"] = list()
		self.elements["O"] = list()

		if chainsToParse == 'all':
			chainsToParseList =  s.get_chains()
		elif chainsToParse == 'firstChainOnly':
			chainsToParseList =  [s.get_chains().next()]
		else:
			chainsToParseList =  []

		for chain in chainsToParseList:
			self.chainEntries[chain.id] = len(self.coords)
			self.chainLen[chain.id] = int()
			for residue in s[0][chain.id]:
				try:
					heteroCode,resNum,_ = residue.get_id()
					if heteroCode == ' ': #if not an HeteroAtom
						self.coords += [residue["CA"].get_coord()]
						self.resID += [(resNum,chain.id)]
						self.residues += [AAVect[AAIndex.index(residue.get_resname()),:]]
						self.sequence += [AAIndex.index(residue.get_resname())]
						self.chainLen[chain.id] += 1
						#print len(self.sequence),(resNum,chain.id),residue.get_resname(),residue["CA"].get_coord()
				except:
					pass
				try:
					heteroCode,resNum,_ = residue.get_id()
					if heteroCode == ' ': #if not an HeteroAtom
						self.elements["N"].append((len(self.coords) - 1,residue["N"].get_coord()))
				except:
					pass
				try:
					heteroCode,resNum,_ = residue.get_id()
					if heteroCode == ' ': #if not an HeteroAtom
						self.elements["F"].append((len(self.coords) - 1,residue["F"].get_coord()))

				except:
					pass
				try:
					heteroCode,resNum,_ = residue.get_id()
					if heteroCode == ' ': #if not an HeteroAtom
						self.elements["S"].append((len(self.coords) - 1,residue["S"].get_coord()))

				except:
					pass
				try:
					heteroCode,resNum,_ = residue.get_id()
					if heteroCode == ' ': #if not an HeteroAtom
						self.elements["O"].append((len(self.coords) - 1,residue["O"].get_coord()))

				except:
					pass



		self.seqLen = len(self.coords)
		self.residues = numpy.array(self.residues)
		self.parsed = True

	# Create Distance Matrix
	def createDM(self,rect = (0,0,0,0)):
		
		if not self.parsed:
			self.parseCoords(self.filename)


		if rect == (0,0,0,0):
			Nstart = 0
			Mstart = 0
			Nend = len(self.coords)
			Mend = len(self.coords)
			self.dmCreated = True
		else:
			Nstart,Nend,Mstart,Mend = rect
			Nend += 1 # rect is defined with all elements included
			Mend += 1 # rect is defined with all elements included


		self.dm = numpy.zeros((len(self.coords),len(self.coords)))
		for i in range(Nstart,Nend):
			for j in range(Mstart,Mend):
				
				self.dm[i,j] = ((self.coords[i][0] - self.coords[j][0])**2 + (self.coords[i][1] - self.coords[j][1])**2 + (self.coords[i][2] - self.coords[j][2])**2)**.5
		


	# Create Distance Matrix
	def createDMStrip(self,width,rect = (0,0,0,0)):
		
		if not self.parsed:
			self.parseCoords(self.filename)

		if rect == (0,0,0,0):
			Nstart = 0
			Mstart = 0
			Nend = len(self.coords)
			Mend = len(self.coords)
		else:
			Nstart,Nend,Mstart,Mend = rect
			Nend += 1 # rect is defined with all elements included
			Mend += 1 # rect is defined with all elements included
			Mend = min(Mend,len(self.coords))
			Nend = min(Nend,len(self.coords))

		self.dm = numpy.zeros((rect[1]-rect[0]+1,rect[3]-rect[2]+1))
		for i in range(Nstart,Nend):
			for j in range(max(0,i-width),min(i+width+1,Mend)):
				self.dm[i,j] = ((self.coords[i][0] - self.coords[j][0])**2 + (self.coords[i][1] - self.coords[j][1])**2 + (self.coords[i][2] - self.coords[j][2])**2)**.5


		return self.dm

	def getHBonds(self):
		if not self.HBondsParsed:
			self.HBonds = list()
			try:
				os.system("dssp -i \"" + self.filename + "\" -o /tmp/tmp.dssp")
				f = open("/tmp/tmp.dssp",'r')
				dsspOutput = f.readlines()

				HBonds = set()
				dsspToPdbID = dict()
				dsspToPdbID[0] = 0
				goParsing = False
				for line in dsspOutput:
					if line.strip() and '#' == line.strip()[0]:
						goParsing = True
					elif goParsing:
						dssp_id = int(line[:5])
						if line[5:11].strip():
							pdb_id = int(line[5:11])
							dsspToPdbID[dssp_id] = dsspToPdbID[dssp_id - 1] + 1
						else:
							dsspToPdbID[dssp_id] = dsspToPdbID[dssp_id - 1]
						bp1 = int(line[25:29])
						bp2 = int(line[29:33])
						if bp1:
							HBonds.add(tuple(sorted([dssp_id,bp1])))
						if bp2:
							HBonds.add(tuple(sorted([dssp_id,bp2])))

				"""
				for Nelement in self.elements["N"]:
					for Oelement in self.elements["O"]:
						if numpy.linalg.norm(Nelement[1] - Oelement[1]) < self.hbondTS:
							vector1 = Nelement[1] - Oelement[1]
							vector2 = self.coords[Oelement[0]] - Oelement[1]
							if math.acos(numpy.dot(vector1,vector2)/(numpy.linalg.norm(vector1) * numpy.linalg.norm(vector2))) * 180 / math.pi < self.hbondAngleTS:
								if abs(Oelement[0] - Nelement[0]) >  1:
									HBonds += [(Oelement[0],Nelement[0])]
				"""
				self.HBonds = [(dsspToPdbID[dsspIdA],dsspToPdbID[dsspIdB]) for (dsspIdA,dsspIdB) in list(HBonds)]
				self.HBondsParsed = True
			except:
				print "OUPS ! Parsing Error  for Hbonds !"
				pass

		return self.HBonds

	def getDM(self):
		if not self.dmCreated:
			self.createDM()
		return self.dm

	# If pdbIdMode is set to true, then we translate back into our internal indexing system
	def getPortionOfDM(self,rectangle,pdbIDMode = False):

		if pdbIDMode:
			(_,chainL),(_,_),(_,chainR),(_,_) = rectangle
			xm,xM,ym,yM = self.getIndexFromResIDRect(rectangle)
		else:
			xm,xM,ym,yM = rectangle

		if not self.dmCreated:
			self.createDM((xm,xM,ym,yM))

		return self.dm[xm:xM+1,ym:yM+1]



	def createCM(self,thsld = DEFAULT_THRESHOLD):
		if not self.parsed: #if not already parsed, parse it
			self.parseCoords(self.filename)
		if not self.dmCreated: #do not need to create two times the distance matrix
			self.createDM()
		self.cm = contactMap(self.dm,thsld)
		return self.cm.cm

	def getCM(self):
		return self.cm.cm

	def loadCM(self,cm):
		self.cm = cm
		self.seqLen = max(max(self.cm.cm[:,0]),max(self.cm.cm[:,1]))+1


	def getCoords(self):
		if not self.parsed:
			self.parseCoords(self.filename)
		return self.coords

	def getSequence(self,start = -1,end = -1):
		if start == -1:
			start = 0
		if end == -1:
			end = len(self.sequence)
		return self.sequence[start:end]




	def getResidues(self):
		return self.residues

	def sort(self,by = 'start'):
		self.parseSS(self.filename)
		if by == 'start':
			self.ss = sorted(self.ss, key = lambda SecStruct : SecStruct.start)
		elif by == 'type':
			self.ss = sorted(self.ss, key = lambda SecStruct : SecStruct.ssType)
		return self


	def getSeqLen(self):
		return self.seqLen

	# get residue number corresponding to the ith atom
	def getResID(self,i):
		return self.resID[i]

	def getResIDRect(self,(xm,xM,ym,yM)):
		return (self.getResID(xm),self.getResID(xM),self.getResID(ym),self.getResID(yM))

	# xm in format (indexxm,chainxm) and so on
	def getIndexFromResIDRect(self,(xm,xM,ym,yM)):
		if xm in self.resID and xM in self.resID and ym in self.resID and yM in self.resID:
			return (self.resID.index(xm),self.resID.index(xM),self.resID.index(ym),self.resID.index(yM))
		else:
			return (float('nan'),float('nan'),float('nan'),float('nan'))

	def __repr__(self):
		return self.ss

	def __str__(self):
		out = ''
		for s in self.ss:
			out += s.__str__() + '\n'
		return out




def sign(x):
	if x>0:
		return 1
	elif x<0:
		return -1
	else:
		return 0



#############################################



# Given a seed (contact) i,j return the smallest rectangle
# including the connected path of (i,j)
def getMinRect((i,j),M,cur = (-1,-1,-1,-1)):
	im = numpy.array(M)
	if cur == (-1,-1,-1,-1):
		cur = (i,i,j,j)

	xm,xM,ym,yM = cur

	if i<xm:
		xm = i
	if j<ym:
		ym = j
	if i>xM:
		xM = i
	if j>yM:
		yM = j

	im[i,j] = -100 # visited

	# W
	if i and im[i-1,j] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j),im,(xm,xM,ym,yM))
	# S
	if j and im[i,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i,j-1),im,(xm,xM,ym,yM))
	# E
	if i < im.shape[0]-1 and im[i+1,j] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j),im,(xm,xM,ym,yM))
	# N
	if j < im.shape[1]-1 and im[i,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i,j+1),im,(xm,xM,ym,yM))

	# SW
	if i and j and im[i-1,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j-1),im,(xm,xM,ym,yM))
	# SE
	if i < im.shape[0]-1 and j and im[i+1,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j-1),im,(xm,xM,ym,yM))
	# NW
	if i and j < im.shape[1]-1 and im[i-1,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j+1),im,(xm,xM,ym,yM))
	# NE
	if i < im.shape[0]-1 and j < im.shape[1]-1 and im[i+1,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j+1),im,(xm,xM,ym,yM))


	return ((xm,xM,ym,yM),im)

# cm is a list of couples in contact
def listRectangles(cm,dmTau):
	if cm:
		# v indicates the visited points by getMinRect of the dist mat by indicating them as '-100'
		(xm,xM,ym,yM),v = getMinRect(cm[0],dmTau)

		# Now we remove all contacts covered by the current connected component
		newCm = list()
		for i,j in cm:
			# If not in the connected component (i.e. not visited by getMinRect)
			if v[i,j] > 0 :
				newCm += [[i,j]]
		return [(xm,xM,ym,yM)] + listRectangles(newCm,dmTau)
				
	else:
		return list()

def getMinMonotonousRect((i,j),M,cur = (-1,-1,-1,-1)):
	im = numpy.array(M)
	if cur == (-1,-1,-1,-1):
		cur = (i,i,j,j)

	xm,xM,ym,yM = cur

	if i<xm:
		xm = i
	if j<ym:
		ym = j
	if i>xM:
		xM = i
	if j>yM:
		yM = j

	im[i,j] = -100 # visited



	if i and im[i-1,j] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j),im,(xm,xM,ym,yM))
	if j and im[i,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i,j-1),im,(xm,xM,ym,yM))
	if i < im.shape[0]-1 and im[i+1,j] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j),im,(xm,xM,ym,yM))
	if j < im.shape[1]-1 and im[i,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i,j+1),im,(xm,xM,ym,yM))


	if i and j and im[i-1,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j-1),im,(xm,xM,ym,yM))

	if i < im.shape[0]-1 and j and im[i+1,j-1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j-1),im,(xm,xM,ym,yM))

	if i and j < im.shape[1]-1 and im[i-1,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i-1,j+1),im,(xm,xM,ym,yM))

	if i < im.shape[0]-1 and j < im.shape[1]-1 and im[i+1,j+1] > 0:
		(xm,xM,ym,yM),im = getMinRect((i+1,j+1),im,(xm,xM,ym,yM))


	return ((xm,xM,ym,yM),im)

# Gives rectangles of CM such that the inner set of 
# contact is either increasing or decreasing
def listMonotonousRectangles(cm,dmTau):
	if cm:
		# v indicates the visited points by getMinRect of the dist mat by indicating them as '-100'
		(xm,xM,ym,yM),v = getMinMonotonousRect(cm[0],dmTau)

		# Now we remove all contacts covered by the current connected component
		newCm = list()
		for i,j in cm:
			# If not in the connected component (i.e. not visited by getMinRect)
			if v[i,j] > 0 :
				newCm += [[i,j]]
		return [(xm,xM,ym,yM)] + listMonotonousRectangles(newCm,dmTau)
				
	else:
		return list()




def isContactInRectangle(contactList,rectangle):
	xm1,xM1,ym1,yM1 = rectangle
	for (a,b) in contactList:
		if xm1<=a and xM1>=a and ym1<=b and yM1>=b:
			return True
		if xm1<=b and xM1>=b and ym1<=a and yM1>=a:
			return True
	return False

# Give the position in the Pad of a contact
def ListContactsRelativesInRectangle(contactList,rectangle,padSize):
	retList = list()
	xm1,xM1,ym1,yM1 = rectangle
	shiftX = (padSize - (xM1 - xm1))/2
	shiftY = (padSize - (yM1 - ym1))/2
	for (a,b) in contactList:
		if xm1<=a and xM1>=a and ym1<=b and yM1>=b:
			retList += [(shiftX + a-xm1,shiftY + b-ym1)]
		if xm1<=b and xM1>=b and ym1<=a and yM1>=a:
			retList += [(shiftX + a-ym1,shiftY + b-xm1)]
	return retList


# First step of Dyn Prog technique to compute the distance between
# two matrices (motifs) using the best alignement. This alignement is 
# constrained to be not gapped farther than Ni Nj along the diagonal
# !! Very Slow !!
def distAlignMotif(dm1,dm2):
	N = dm1.shape[0]
	scores = - numpy.ones((N,N,N,N))
	path = list()
	maxScore = 0

	gapScore = 10
	Ni = 3
	Nj = 3


	
	for i in range(N):
		for j in range(N):
			for k in range(max(1,i-Ni),min(N,i+Ni)):
				for l in range(max(1,j-Nj),min(N,j+Nj)):
					tmpScore = -1
					tmpDirection = (0,0,0)
					# sx indicates the shift value of var x
					for si in [-1,0]:
						for sj in [-1,0]:
							for sk in [-1,0]:
								for sl in [-1,0]:
									if si+sj+sk+sl<0: # in case something is shifted
										# we take the min value of the predecesors
										if  scores[max(0,i+si),max(0,j+sj),k+sk,l+sl]>0 and (tmpScore < 0 or tmpScore > scores[max(0,i+si),max(0,j+sj),k+sk,l+sl]):
											tmpScore  = scores[max(0,i+si),max(0,j+sj),k+sk,l+sl]
											tmpDirection = (max(0,si),max(0,sj),sk,sl)

					scores[i,j,k,l] = tmpScore + (dm1[i,j] - dm2[k,l])**2 # + gapScore*(4-sum(tmpDirection))
					maxScore  = max(maxScore,scores[i,j,k,l])


	return numpy.sqrt(maxScore)


def modulesOfArrayElements2d(m):
	modules = numpy.zeros(m.shape)
	for i in range(m.shape[0]):
		for j in range(m.shape[1]):
			modules[i,j] = abs(m[i,j])
	return modules



def getEnergySpectra2d(m):
	return modulesOfArrayElements2d(fft2(m))



def loadSequencePair(fileName):
	f = open(fileName, 'rb')
	buf = f.readlines()
	seqs = [eval(buf[0]),eval(buf[1])]
	indices = eval(buf[2])
	chains = (buf[3][0],buf[4][0])
	f.close()
	return seqs,indices,chains

def writeSequencePair(fileName,(a,b),(xm,xM,ym,yM), (chainA,chainB)):
	f = open(fileName, 'wb')
	f.write(str(a) + '\n')
	f.write(str(b) + '\n')
	f.write(str([xm,xM,ym,yM]) + '\n')
	f.write(str(chainA) + '\n')
	f.write(str(chainB) + '\n')
	f.close()

def WriteFASTAFile(fileName,pool,motifNames):
	FastaSeq = '\n'.join(['>' + motifNames[i] + '\n' + str().join([AANick[AAIndex[x]] for x in pool[i]]) for i in range(len(pool))])
	f = open(fileName, 'w')
	f.write(FastaSeq)
	f.close()

def getResiduesInSequencePair(p,i,s):
	a = p[s]
	if abs(i-padSize/2)>=len(a)/2:
		return -1
	else:
		return a[i+len(a)/2-padSize/2]

########################
# Utils
########################
def present(l,item):
	try:
		i = l.index(item)
	except ValueError:
		i = -1 # no match
	return True if i>=0 else False


def mult((a,b)):
	return a*b


def padMyMatrix(m,p):
	l,h = m.shape
	paddedMatrix = numpy.zeros((p,p))
	paddedMatrix[(p-l)/2:(p+l)/2,(p-h)/2:(p+h)/2] = m
	return paddedMatrix















