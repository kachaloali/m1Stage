%* fullclassify4.lp
 	Assignment of the most likely class to some test sequences from their paloma blocks composition
	given the block's composition of learning sequences with known class
	Author: Jacques Nicolas

 One consider only blocks that  appear both in a test sequence and in another sequence 
 The formal concepts BlocksXSequences for each class (including the class unknown of test sequences)
 are produced using program blockconcept.lp

 The goal is to propose a classification for each test sequence
	
%Data:
 sequence(Ident seq, Class).
 dubious(Learning Sequence Seq, Declared dubious Class of Seq, Most likely Class of Seq).
 concept(N° concept, Subset type -extent for Blocks or intent for Sequences, Subset element)


%Results:

 - bestclassified(S,C,NB,NS):  test sequence S has a unique most specific classification in class C 
 and it is suported by NB blocks and NS sequences having these blocks in common with S
 - ambiguous(S,C1,NB1,NS1): test sequence S has at least two most specific classifications,  
 one being class C that is suported by NB blocks and NS sequences having these blocks in common with S
 - new_family(I): There exists a new family corresponding to concept I,
 that is, made of a subset of sequences that have a characteristic set of common blocks described in I
 - marker_new-family(S): the test sequence S that is both classifed and member of a new family
 - support_new_family(I,NB): a new family is described by concept I that contains NB blocks
 - unclassified(S): test sequence S can neither be classified nor be associated to the existence of a new family

%Call
 cd C:\Users\jnicolas\Documents\potsdam\ASP\paloma
  ..\bin\clingopyth52.exe --stats resblockclassconcept.lp fullclassify5.lp
  ..\bin\clingopyth52.exe --stats Notes_S3_t1m1M10cAlignX-concepts.txt  fullclassify5.lp
*%

%Class of unlabeled instances to be assigned with a class label
#const unknown=black.

%Maximum allowed number of new families
#const maxnew=20.
%Maximum number of classes for a new family chimera
#const maxchimeraclass=2.

param(list(unknown(unknown), maxnew(maxnew), maxchimeraclass(maxchimeraclass))).

%Test sequences in each concept
test(Seq):- sequence(Seq,unknown).
%test(I,Seq):- concept(I,intent,Seq), test(Seq).

%Learning sequences in each concept
learning(I,Seq):- concept(I,intent,Seq), sequence(Seq,C), C!=unknown.
learning(Seq):- learning(_,Seq).

%Concept numbers
num(I):- concept(I,_,_).

%Classes present in concepts
class(C):- concept(_,intent,C), sequence(_,C).
classconcept(I,C):- concept(I,intent,C), class(C).


%Block assciated to some class in a concept
blockclass(B,C):- concept(I,extent,B), classconcept(I,C).


%Number of blocks for each concept
block_support(I,N):- N={concept(I,extent,B)}, num(I).

%Number of learning sequences for each concept
learning_seq_support(I,N):-   learning(I,_), N={concept(I,intent,S):learning(S)}, num(I).

%Blocks of a concept with a unique class 
characteristic_block(I,Class,B):- concept(I,extent,B), blockclass(B,Class), not blockclass(B,C):class(C), C!=Class.
characteristic_concept(I,Class):- characteristic_block(I,Class,_).

%Number of characteristic blocks in each concept 
characteristic_block_support(I,N):- characteristic_concept(I,_), N={characteristic_block(I,Class,B)}, num(I).



%Concepts with the same non empty subset of associated classes
o_same_classes(I,J):-  num(I); num(J); I<J; classconcept(I,D); classconcept(J,D); classconcept(I,C):classconcept(J,C); classconcept(J,C):classconcept(I,C).
same_classes(I,J):- o_same_classes(I,J).
same_classes(I,J):- o_same_classes(J,I).
same_classes_or_equal(I,J):-  same_classes(I,J).
same_classes_or_equal(I,I):-  num(I).


%Concepts with the same subset of associated test sequences
%o_same_test(I,J):- num(I); num(J); I<J; concept(I,intent,S):concept(J,intent,S),test(S); concept(J,intent,S):concept(I,intent,S),test(S).
%same_test(I,J):- o_same_test(I,J).
%same_test(I,J):- o_same_test(J,I).


% The concept I  contains a strict subset of blocks of  concept J, 
% I and J shares a common test sequence S 
blocksincluded(I,J,S):- concept(I,intent,S), concept(J,intent,S), test(S), classconcept(I,C), classconcept(J,C), I!=J,
	concept(J,extent,B):concept(I,extent,B).
	
% Concepts with a maximal number of blocks (subset maximization) including both a test sequence S and a class C 
maximal_concept(S,I,C):- concept(I,intent,S),test(S), classconcept(I,C), not blocksincluded(I,J,S): blocksincluded(I,J,S).


% Maximal concept with a largest block support of size NB including both a test sequence S and a class C 
top_concept(S,I,C,NB):- maximal_concept(S,I,C), block_support(I,NB), 
	not top_concept(S,J,U,X):maximal_concept(S,J,U),block_support(J,X),X>NB. 
top_concept(S,I,C):-top_concept(S,I,C,_).

% Maximal concept with a largest learning sequence support  including both a test sequence S and a class C 
% containing a given test sequence S, a given class C and the same classes than the corresponding top concept
bottom_concept(S,J,C,NS):- top_concept(S,I,C),  same_classes_or_equal(I,J),  learning_seq_support(J,NS), 
	not bottom_concept(S,K,C,X):same_classes(J,K),learning_seq_support(K,X),X>NS. 

%Characteristic concept with a largest set of blocks containing a given test sequence S and a class C. 
%N is the number of characteristic blocks
possible_characteristic_top_concept(S,I,C,N):- top_concept(S,I,C), characteristic_concept(I,C), 
	characteristic_block_support(I,N).

characteristic_top_concept(S,I,C,NB):- 	possible_characteristic_top_concept(S,I,C,N); top_concept(S,I,C,NB);
	not possible_characteristic_top_concept(S,J,U,X):maximal_concept(S,J,U),characteristic_block_support(J,X),X>N. 

characteristic_bottom_concept(S,C,NS):- characteristic_top_concept(S,I,C,NB); same_classes(I,J);  learning_seq_support(J,NS); 
	not bottom_concept(S,K,C,X):same_classes(J,K),characteristic_concept(K,C),learning_seq_support(K,X),X>NS; 
	characteristic_concept(I,C).



%The test sequence Seq is specific of class C in concept I and has support of NB blocks and NS sequences for that.
specific_class(Seq,I,C,NB,NS):- characteristic_top_concept(Seq,I,C,NB), characteristic_bottom_concept(S,C,NS).
specific_class(Seq,C):- specific_class(Seq,_,C,_,_).

%The test sequence Seq is specific of class Class in all concepts and has support of NB blocks and NS sequences for that.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bestclassified(Seq,Class,I,NB,NS):- specific_class(Seq,I,Class,NB,NS), not specific_class(Seq,C):class(C),C!=Class.
bestclassified(Seq):- bestclassified(Seq,_,_,_,_).



%Class compatible with a test sequence Seq and has support of NB blocks and NS sequences for that.
max_compatible_class(Seq,C,I,NB):- top_concept(Seq,I,C,NB), not bestclassified(Seq).
max_compatible_class(Seq,C):- max_compatible_class(Seq,C,_,_).

%A test sequence is classified if it is compatible with only one class and not already bestclassified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classified(Seq,C,I,NB,NS):- max_compatible_class(Seq,C,I,NB),  bottom_concept(Seq,_,C,NS), not max_compatible_class(Seq,C2):class(C2),C2!=C.
classified(S):- classified(S,_,_,_,_).


%The test sequence Seq is specific of at least two classes.
%For class Class in concept I, it has support of NB blocks and NS sequences for that.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bestambiguous(Seq,Class,I,NB,NS):- specific_class(Seq, I, Class, NB,NS), not bestclassified(Seq), not classified(Seq).
bestambiguous(Seq,Class):- bestambiguous(Seq,Class,I,NB,NS).
bestambiguous(Seq):- bestambiguous(Seq,_).


%A test sequence is ambiguous if it can be predicted in two different classes and is not already bestambiguous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ambiguous(Seq,C1,I1,NB1,NS1):- max_compatible_class(Seq,C1,I1,NB1),  bottom_concept(Seq,_,C,NS1), 
	max_compatible_class(Seq,C2,I2,NB2), C1!=C2, not bestambiguous(Seq).
ambiguous(Seq,C1):- ambiguous(Seq,C1,I1,NB1,NS1).
ambiguous(S):- ambiguous(S,_).

%The test sequence Seq is specific of at least two new classes. Not possible here since new classes are not defined
%bestnewcomposition(Seq,I):- bestambiguous(Seq,Class,I,NB,NS), learning_seq_support(I,0).

%newcomposition(Seq,I):- ambiguous(Seq,Class,I,NB,NS), learning_seq_support(I,0).
%newcomposition(Seq):- newcomposition(Seq,I).


%
%Choice of concepts corresponding to new families
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Concepts and Blocks associated to new families (with an unknown class, not represented)
new_family_block(I,B):- concept(I,extent,B), not concept(I,intent,C):class(C).
new_family_block(B):- new_family_block(_,B).
new_family_num(I):- new_family_block(I,_).

%Blocks that are exclusive of test sequences for all concepts
new_test_block(B):- new_family_block(B), new_family_block(I,B):concept(I,extent,B).

unknown(S):- test(S), not bestclassified(S), not classified(S), not bestambiguous(S), not ambiguous(S).

%A sequence of a new family is a test sequence of a new family concept
new_family_seq(S):- new_family_num(I), concept(I,intent,S), unknown(S).

%New family with a never observed block B 
new_family_new_block(I,B):- new_test_block(B), concept(I,extent,B).

%%%%
%Choice of a set of new families covering all unknown sequences
%%%%
{new_family(I):new_family_num(I)} maxnew.


%Characteristic blocks and test sequences that belong to chosen new family concepts
blocknew(I,B):- new_family(I), new_family_block(I,B).

membernew(S,I):- new_family(I), concept(I,intent,S), test(S).
membernew(S):- membernew(S,_).

%All new family sequences have to be covered by a new family concept
:-new_family_seq(S), not membernew(S). 


%Sequences in the chosen new family concepts that can be classified without introducing a new class
marker_new(S):- membernew(S,I), not new_family_seq(S).

%At least one new family sequence by concept associated to a new family
:-new_family(I); marker_new(S):concept(I,intent,S).


%ambiguous family with its blocks appearing independently in several classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chimera(S):- 2{bestambiguous(S,C)}maxchimeraclass, test(S).
chimera(S):- 2{ambiguous(S,C)}maxchimeraclass, test(S).  %, not newcomposition(S).
chimera(S,C):- chimera(S), ambiguous(S,C).

%Block support of each new family (by definition, no seq support)
support_new_family(I,NB):- block_support(I,NB), new_family(I).

unclassified(S):-  unknown(S), not membernew(S).


%Minimize the number of already classified sequences in new families
#minimize{1@5,S:marker_new(S)}.
%Minimize the number of concepts associated to new families
#minimize{1@3,I:new_family(I)}.
%Maximize the number of characteristic blocks associated to new families
#maximize{1@1,I,B:blocknew(I,B)}.



#show param/1.
%#show dubious/3.
%#show classconcept/2.
%#show maximal_concept/3.
%#show bottom_concept(sYidA,I,C,NS).
%#show new_test_block/1.
%#show characteristic_top_concept(sYaeD,I,C,NB).
#show bestclassified/5. #show bestclassified/1.
#show bestambiguous/5. #show bestambiguous/1.
#show classified/5.
#show classified/1.
#show ambiguous/5. #show ambiguous/1.
%#show bestnewcomposition/2.
%#show newcomposition/2.
#show marker_new/1.
#show new_family/1.
#show membernew/2.
#show blocknew/2.
%#show new_family_new_block/2.
#show chimera/2.
#show new_family_seq/1.
#show support_new_family/2.
#show unclassified/1.
#show characteristic_block/3.
