%* blockclassconcept2.lp
	Find  formal concepts made of blocks and sequences 
	for each class having a common block with a test sequence
	Run with the clasp option -n0 to have all the formals concepts
	Author: Jacques Nicolas

%Call
 cd C:\Users\jnicolas\Documents\potsdam\ASP\paloma
 ..\bin\clingopyth52.exe blockclassconcept3.lp Notes_S3_t5m1M10CAlignX-logic-out.txt --heu=Vsids --stats -n 0 

%Data:
 sequence(Ident seq, Class).
 block(N° block, Ident seq).

*%
 
#script (python)

import clingo
fichier = open("resblockclassconcept.lp", "w")

def on_model(conc, m):
   atoms= m.symbols(atoms=True)
   for atom in atoms:
      if conc==1:
         if   atom.name == "sequence" and len(atom.arguments) == 2:
            fichier.write( "sequence("+str(atom.arguments[0])+", "+ str(atom.arguments[1])+").\n")
         elif atom.name == "dubious" and len(atom.arguments) == 3 :
            fichier.write( "dubious("+str(atom.arguments[0])+", "+ str(atom.arguments[1])+", "+ str(atom.arguments[2])+").\n")
      if   atom.name == "min_element_intent" and len(atom.arguments) == 1:
         fichier.write("min_element_intent("+str(conc)+", "+str(atom.arguments[0])+").\n")
      elif atom.name == "concept" and len(atom.arguments) == 2 :
         fichier.write("concept("+str(conc)+", "+str(atom.arguments[0])+", "+ str(atom.arguments[1])+").\n")
 


def main(prg):
   prg.ground([("base", [])])
   for ind, answer in enumerate(prg.solve(yield_=True, async=True), start=1):
      on_model(ind, answer)
#      fichier.write("\n")
   fichier.close()
#end.



%Class of unlabelled instances to be assigned with a class label
#const unknown=black.
% Minimum support for concepts
#const minintent=1.
#const minextent=1.


%-----------------------------------------------------
% General definitions
%-----------------------------------------------------

param(list(unknown(unknown))).

%List of all classes present in data
	classu(Class):- sequence(_, Class).
%List of all classes except the unknown class (learning classes)
	class(Class):- classu(Class), Class!=unknown.
%Number of classes present in data
	nbclasses(N):-N={classu(Class)}.
%List of Classes of a given block
	class(Block,Class):- block(Block, Seq), sequence(Seq,Class).

%Test sequences -to be to be assigned with a class label-
	test(Seq):- sequence(Seq, unknown).
%Number of test sequences
	nbtest(N):-N={test(Seq)}.
%Sequences of the learning set
	learningset(Seq):- sequence(Seq, C), class(C).
%Number of learning sequences
	nblearningseq(N):-N={learningset(Seq)}.

	
%---------------------------------------------------------------------------------------
% Analysis of the quality of class annotations
% A learning sequence is considered dubious if it is classified in some class Class1
% but the number of "interesting" blocks for some class Class2 differing from Class1
% is greater than for the cumulated sum of "interesting" blocks over all other classes.
%---------------------------------------------------------------------------------------

%Interesting blocks of a sequence are either blocks of a class (present in all its sequences) 
% or blocks characteristic of a class (present uniquely in this class apart possibly the sequence)
% interesting_block(Sequence Seq, Block of Seq, Class)
	interesting_block(Seq,Block,Class):- learningset(Seq), block(Block, Seq), classu(Class),  block(Block, S):sequence(S,Class).
	interesting_block(Seq,Block,Class):- learningset(Seq), block(Block, Seq), classu(Class),  sequence(S,Class):block(Block, S),S!=Seq.

% dubious (Strangely Classified Learning Sequence Seq, Class of Seq, Most likely Class of Seq).
	dubious(Seq,Class1,Class2):- learningset(Seq), sequence(Seq, Class1), class(Class2), Class2!=Class1,
		M={interesting_block(Seq,_,C): class(C), C!=Class2}, M+1{interesting_block(Seq,_,Class2)}.
	dubious(Seq):- dubious(Seq,Class1,Class2).

%an admissible sequence is a test sequence or non dubious learning sequence
	admissible_seq(Seq):-  learningset(Seq), not dubious(Seq).
	admissible_seq(Seq):-  test(Seq).

	
%-----------------------------------------------------
% Build the binary relation for Formal concept analysis
%-----------------------------------------------------

%Common block between a test sequence Seq and an admissible  sequence S
	validblock(Block):- block(Block, Seq), test(Seq), block(Block, S), admissible_seq(S), S!=Seq.

%Binary relation including the Classes: a class contains a block if at least one of its sequence contains it
	rel(Block, Seq):- block(Block, Seq), validblock(Block).
	rel(Block, Class):- rel(Block, Seq), learningset(Seq), sequence(Seq,Class).

	extent(Extent):- rel(Extent, _).
	intent(Intent):- rel(_, Intent).
	

---------------------------
% Find the formal concepts
%--------------------------

% Choose a non empty intent
	minintent{ concept(intent,X):intent(X) }.

% Minimum element of the intent subset of the concept
	min_element_intent(Min):-  concept(intent,Min),  not concept(intent,X):intent(X),X<Min.
	
% Choose a non empty subset of extent that are related to the chosen intent
	minextent{ concept(extent,Y):rel(Y,M) } :-  min_element_intent(M).

% Extent that is not connected to some intent in the concept
	imp_extent(X) :-  extent(X), concept(intent,Y), not rel(X,Y).
% Intent that is not connected to some extent in the concept
	imp_intent(Y) :-  intent(Y), concept(extent,X),  not rel(X,Y).

% Consistency: no impossible extent or intent belonging to the concept
	:- imp_extent(X), concept(extent,X).
	:- imp_intent(Y), concept(intent,Y).

% Maximality: all possible extent or intent have to be added to the concept
	:- not imp_extent(X), not concept(extent,X), extent(X).
	:- not imp_intent(Y), not concept(intent,Y), intent(Y).


#show nbtest/1.
#show nblearningseq/1.
