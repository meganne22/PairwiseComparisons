(* ::Package:: *)

(* Pairwise Comparison (Ranking) Package for Mathematica(TM) *)
(* (C) 2014 Konrad Kulakowski *)
(* e-mail: konrad@kulakowski.org *)

(* version 0.2 *)

BeginPackage["PairwiseComparisons`"]

(**********************************************************************************************************)
(* Eigenvalue Rank Methods *)
(**********************************************************************************************************)

PrincipalEigenValue::usage = "PrincipalEigenValue[M] returns the largest principal eigenvalue of M";
PrincipalEigenVector::usage = "PrincipalEigenVector[M] returns the eigenvector of M corresponding to its principal eigenvalue ";
SaatyIdx::usage = "SaatyIdx[M] returns the value of the Saaty Inconsistency Index computed for the matrix M";
EigenvalueRank::usage = "EigenvalueRank[M] returns the principal eigen vector of M rescalled in a way that the sum of its entries is 1";

(**********************************************************************************************************)
(* Geometric Mean Rank Methods *)
(**********************************************************************************************************)

GeometricRank::usage = "GeometricRank[M] returns rank list given as geometric means of rows of the matrix M";
GeometricRescaledRank::usage = "GeometricRescaledRank[M] returns a geometric mean rank list of M rescalled in a way that the sum of its entries is 1";

(**********************************************************************************************************)
(* HRE  (Heuristic Rating Estimation) - rankings with the reference values *)  
(**********************************************************************************************************)

HREConstantTermVector::usage = "HREConstantTermVectorEntry[M,v] returns the constant term vector for the HRE method. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREMatrix::usage = "A = HREMatrix[M,v] returns a matrix for the HRE method. Togeter with b=HREConstantTermVector[M, v] it forms the linear equation system Au=b where u is the partial ranking vector. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREPartialRank::usage = "HREPartialRank[M,v] returns values for the unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREFullRank::usage = "HREFullRank[M,v] returns values for both: known and unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HRERescaledRank::usage = "HRERescaledRank[M,v] returns the rescaled HREFullRank list, rescalled in a way that the entries sum up to 1"

(**********************************************************************************************************)
(* HRE Geometric (Heuristic Rating Estimation using Geometric Mean) - rankings with the reference values *)  
(**********************************************************************************************************)

HREGeomConstantTermVector::usage = "HREGeomConstantTermVectorEntry[M,v] returns the constant term vector for the HRE method. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREGeomMatrix::usage = "A = HREGeomMatrix[M,v] returns a matrix for the HRE Geometric method. Togeter with b=HREGeomConstantTermVector[M, v] it forms the linear equation system Au=b where u is the intermediate partial ranking vector. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREGeomIntermediateRank::usage = "HREGeomIntermediateRank[M,v] returns values for the unknown alternatives before they are raised to the power 10. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREGeomPartialRank::usage = "HREGeomPartialRank[M,v] returns values for the unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREGeomFullRank::usage = "HREGeomFullRank[M,v] returns values for both: known and unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"
HREGeomRescaledRank::usage = "HREGeomRescaledRank[M,v] returns the rescaled HREGeomFullRank list, rescalled in a way that the entries sum up to 1"

(**********************************************************************************************************)
(* Koczkodaj Inconsistency Methods *)
(**********************************************************************************************************)

KoczkodajTriadIdx::usage = "KoczkodajTriadIdx[{a,b,c}] returns the Koczkodaj triad inconsistency i.e. Min[Abs[1-c/ab],Abs[1-ab/c]]";
KoczkodajIdx::usage = "KoczkodajIdx[M] returns the value of the Koczkodaj inconsistency computed for the matrix M";
KoczkodajTheWorstTriad::usage = "KoczkodajTheWorstTriad[M] returns the worst triad in M according to Koczkodaj inconsistency criterion";
KoczkodajTheWorstTriads::usage = "KoczkodajTheWorstTriads[M, n] returns n worst triads in M according to Koczkodaj inconsistency criterion";
KoczkodajConsistentTriad::usage = "KoczkodajConsistentTriad[{a,b,c}] returns the closest consistent triad";
KoczkodajImproveMatrixStep::usage = "KoczkodajImproveMatrixStep[M] returns improved (more consistent) matrix in which the most inconsistent triad is replaced by the closest consistent one";

(**********************************************************************************************************)
(* Condition Of Order Preservation by Bana e Costa and Vansnick *)
(**********************************************************************************************************)

COP1ViolationList::usage = "COP1ViolationList[M,r] returns the list of indices that violate the first Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick). Here, M is the PC Matrix, whilst r is the ranking result list."

COP1Check::usage = "COP2ViolationList[M,r] it is true when the first Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick) is violated. Here, M is the PC Matrix, whilst r is the ranking result list."

COP2ViolationList::usage = "COP2ViolationList[M,r] returns the list of indices that violate the second Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick). Here, M is the PC Matrix, whilst r is the ranking result list."

COP2Check::usage = "COP2ViolationList[M,r] it is true when the second Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick) is violated. Here, M is the PC Matrix, whilst r is the ranking result list."

(**********************************************************************************************************)
(* Ranking Discrepancy  *)
(**********************************************************************************************************)

ErrorMatrix::usage = "ErrorMatrix[M, mju] computes matrix containing entries in form e_ij = m_ji*mju(c_i)/mju(c_j). When the matrix is consistent every e_ij equals 1"

LocalDiscrepancyMatrix::usage = "LocalDiscrepancyMatrix[M, mju] computes matrix with entries d_ij = max{e_ij - 1, 1/e_ij - 1}"

GlobalDiscrepancy::usage = "GlobalDiscrepancy[M, mju] returns maximal entry of  LocalDiscrepancyMatrix[M, mju]"

(**********************************************************************************************************)
(* Auxiliary functions *)
(**********************************************************************************************************)

RecreatePCMatrix::usage = "RecreatePCMatrix[M] recreates reciprocal matrix on the basis of upper-triagle of M";
GetMatrixEntry::usage = "GetMatrixEntry[i ,j, M] - the same as M[i,j]"
DeleteRows::usage = "TODO"
DeleteColumns::usage = "TODO"
DeleteRowsAndColumns::usage = "TODO"
SetDiagonal::usage = "TODO"
RankToVector::usage = "TODO"
RankOrder::usage = "RankOrder[rankList] - prints order of elements on the rank list"

Begin["`Private`"]

(***************** Eigenvalue based rankings - Implementation *****************) 

PrincipalEigenValue[matrix_] := Max[Abs[N[Eigenvalues[matrix]]]]

PrincipalEigenVector[matrix_] := 
 Extract[N[Eigenvectors[matrix]], 
  First[First[
    Position[Abs[N[Eigenvalues[matrix]]], 
     PrincipalEigenValue[matrix]]]]]
 
SaatyIdx[matrix_] := With[{n = First[Dimensions[matrix]], alpha = PrincipalEigenValue[matrix]}, (alpha - n)/(n - 1)]

EigenvalueRank[matrix_] := With[{pev = PrincipalEigenVector[matrix]}, (#/Apply[Plus, pev]) & /@ pev]

(***************** Geometric mean rankings - Implementation *****************)

GeometricRank[matrix_] := GeometricMean[#] & /@ matrix

(* deprecated *)  
NGeometricRank[matrix_] := N[GeometricMean[#] & /@ matrix]

GeometricRescaledRank[matrix_] := With[{gr = GeometricRank[matrix]}, (#/Apply[Plus, gr]) & /@ gr ]

(* deprecated *)  
NGeometricRescaledRank[matrix_] := N[GeometricRescaledRank[matrix]]
 
(***************** HRE - Rankings with the reference values - Implementation *****************)

GetMatrixEntry[i_, j_, matrix_] := Extract[Extract[matrix, i], j]

auxMakeRecip[a_, b_, matrix_] :=
 If[Extract[b, 1] < Extract[b, 2], a,
  If[Extract[b, 1] == Extract[b, 2], 1, 
   1/GetMatrixEntry[Extract[b, 2], Extract[b, 1], matrix]]]
  
RecreatePCMatrix[upperTriangleMatrix_] := 
 MapIndexed[auxMakeRecip[#1, #2, upperTriangleMatrix] &, 
  upperTriangleMatrix, {2}]

zeroQ[x_] := (ToString[Head[x]] == "Integer" ||  ToString[Head[x]] == "Real") && x == 0 
zeroIndices[knownVectorList_] :=  Map[First[Last[#]] &, Select[MapIndexed[List, knownVectorList], zeroQ[First[#]] &]] 
nonZeroValues[knownVectorList_] := Select[knownVectorList, Not[zeroQ[#]] &]
nonZeroIndices[knownVectorList_] := Map[First[Last[#]] &,  Select[MapIndexed[List, knownVectorList], Not[zeroQ[First[#]]] &]]
 
HREConstantTermVectorEntry[matrix_, knownVectorList_, j_] := 
 Apply[Plus, 
  MapThread[(#1 * #2 * (1/(First[Dimensions[matrix]] - 1))) &, 
   List[nonZeroValues[knownVectorList], 
    Map[GetMatrixEntry[j, #, matrix] &, 
     nonZeroIndices[knownVectorList]]]]]

HREConstantTermVector[matrix_, knownVector_] := 
 Map[List, 
  Map[HREConstantTermVectorEntry[matrix, Flatten[knownVector], #] &, 
   zeroIndices[Flatten[knownVector]]]]

DeleteRows[matrix_, listOfRows_] :=  Delete[matrix, Map[List, listOfRows]]

DeleteColumns[matrix_, listOfColumns_] := Map[Delete[#, Map[List, listOfColumns]] &, matrix]

DeleteRowsAndColumns[matrix_, listOfIndicesToRemove_] :=  
		   DeleteColumns[DeleteRows[matrix, listOfIndicesToRemove],   listOfIndicesToRemove]

SetDiagonal[matrix_, diagVal_] := MapIndexed[If[First[#2] == Last[#2], diagVal, #1] &, matrix, {2}]
 
HREMatrix[matrix_, knownVector_] := 
 SetDiagonal[(-1/(First[Dimensions[matrix]] - 1))*
   DeleteRowsAndColumns[matrix, nonZeroIndices[Flatten[knownVector]]], 1]

HREPartialRank[matrix_, knownVector_] := 
 LinearSolve[HREMatrix[matrix, knownVector], 
  HREConstantTermVector[matrix, knownVector]]

indexedUnknownValues[hrePartialResult_, knownVector_] := 
 MapThread[
  List, {Flatten[hrePartialResult], 
   zeroIndices[Flatten[knownVector]]}]

indexedKnownValues[knownVector_] := 
 MapThread[
  List, {nonZeroValues[Flatten[knownVector]], 
   nonZeroIndices[Flatten[knownVector]]}]

indexedValues[hrePartialResult_, knownVector_] := 
 Sort[Join[indexedUnknownValues[hrePartialResult, knownVector], 
   indexedKnownValues[knownVector]], #1[[2]] < #2[[2]] &]

HREFullRank[matrix_, knownVector_] := 
 Map[First[#] &, 
  indexedValues[HREPartialRank[matrix, knownVector], knownVector]]

HRERescaledRank[matrix_, knownVector_] :=
  With[{fullRank = HREFullRank[matrix, knownVector]}, (#/Apply[Plus, fullRank]) & /@ fullRank]

(***************** HRE Geom - Rankings with the reference values - Implementation *****************)

HREGeomMatrix[matrix_, knownVector_]:= With[{dim1 = First[Dimensions[matrix]], 
	dim2 = Length[Select[Flatten[knownVector], zeroQ[#] & ]]},
		Map[If[#==0,-1,#]&,(dim1-1)*IdentityMatrix[dim2], {2}]] 

HREGeomConstantTermVectorEntry[matrix_, knownVector_, j_]:=With[{jthRow = Drop[Extract[matrix,{j}], {j}], 
	knownVectEntr = Select[Flatten[knownVector], Not[zeroQ[#]]  & ]},
	Log10[Times @@ Join[jthRow, knownVectEntr]]]

HREGeomConstantTermVector[matrix_, knownVector_] := 
 Map[List, Map[HREGeomConstantTermVectorEntry[matrix, Flatten[knownVector], #] &, 
   zeroIndices[Flatten[knownVector]]]]

HREGeomIntermediateRank[matrix_, knownVector_]:=LinearSolve[HREGeomMatrix[matrix, knownVector],
	HREGeomConstantTermVector[matrix, knownVector]]

HREGeomPartialRank[matrix_, knownVector_]:= Map[10^# &, HREGeomIntermediateRank[matrix, knownVector]]

HREGeomFullRank[matrix_, knownVector_] := 
 Map[First[#] &, indexedValues[HREGeomPartialRank[matrix, knownVector], knownVector]]

HREGeomRescaledRank[matrix_, knownVector_] := With[{fullRank = HREGeomFullRank[matrix, knownVector]},
	Map[#/(Plus @@ fullRank) &, fullRank]]

(***************** Koczkodaj inconsistency methods - Implementation *****************)

KoczkodajTriadIdx[triad_] := 
 Min[Abs[1 - (Extract[triad, {1}]*Extract[triad, {2}])/
     Extract[triad, {3}]], 
  Abs[1 - Extract[
      triad, {3}]/(Extract[triad, {1}]*Extract[triad, {2}])]]    

auxUniqueTuples[n_] := 
 Select[Tuples[Range[n], 3], 
  Extract[#, 1] != Extract[#, 2] &&   
    Extract[#, 1] != Extract[#, 3] &&  
    Extract[#, 2] != Extract[#, 3] &] 
 
auxMakeATriad[matrix_, tuple_] := 
 List[Extract[Extract[matrix , Extract[tuple, 1]], 
   Extract[tuple, 2]], 
  Extract[Extract[matrix , Extract[tuple, 2]], Extract[tuple, 3]], 
  Extract[Extract[matrix , Extract[tuple, 1]], Extract[tuple, 3]]] 
 
(*
auxUniqueTriads[matrix_] := 
 auxMakeATriad[matrix, #] & /@ 
  auxUniqueTuples[First[Dimensions[matrix]]] 

KoczkodajIdx[matrix_] := 
 Max[KoczkodajTriadIdx[#] & /@ auxUniqueTriads[matrix]]  
*)

auxUniqueTriadsAndTuples[matrix_] := 
 List[#, auxMakeATriad[matrix, #], 
    KoczkodajTriadIdx[auxMakeATriad[matrix, #]]] & /@ 
  auxUniqueTuples[First[Dimensions[matrix]]]

KoczkodajTheWorstTriad[matrix_] := 
 First[Sort[auxUniqueTriadsAndTuples[matrix], 
   Last[#1] > Last[#2] & ] ] 

KoczkodajTheWorstTriads[matrix_, n_] := 
 Take[Sort[auxUniqueTriadsAndTuples[matrix], Last[#1] > Last[#2] & ] ,
   n]

KoczkodajIdx[matrix_] := Last[KoczkodajTheWorstTriad[matrix]]

KoczkodajConsistentTriad[triad_] := 
 With[{q = First[triad], r = Extract[triad, {2}], s = Last[triad]},
  List[q^(2/3) r^(-1/3) s^(1/3), q^(-1/3) r^(2/3) s^(1/3), 
   q^(1/3) r^(1/3) s^(2/3)]]

KoczkodajImproveMatrixStep[matrix_] := 
 With[{worstTriad = KoczkodajTheWorstTriad[matrix]}, 
  With[{improvedTriad = KoczkodajConsistentTriad[Extract[worstTriad, {2}]], 
    mik = Drop[First[worstTriad], {3}], 
    mkj = Drop[First[worstTriad], {1}], 
    mij = Drop[First[worstTriad], {2}]}, 
   With[{mikval = First[improvedTriad], 
     mkjval = Extract[improvedTriad, {2}], 
     mijval = Last[improvedTriad]}, 
    ReplacePart[
     ReplacePart[
      ReplacePart[
       ReplacePart[
        ReplacePart[ReplacePart[matrix, Reverse[mij] -> 1/mijval], 
         mij -> mijval], Reverse[mkj] -> 1/mkjval], mkj -> mkjval], 
      Reverse[mik] -> 1/mikval], mik -> mikval]]]]

(***************** Condition of order preservation methods - Implementation *****************)

COP1ViolationList[matrix_, resultList_] := 
 Select[Flatten[
   MapIndexed[ 
    List[((#1 <= 
           1) && (Extract[resultList, First[#2]]/
            Extract[resultList, Last[#2]] <= 1)) || ((#1 >= 
           1) && (Extract[resultList, First[#2]]/
            Extract[resultList, Last[#2]] >= 1)), #2] &, matrix, {2}],
    1], First[#] == False &]

COP1Check[matrix_, resultList_] := 
 Length[COP1ViolationList[matrix, resultList]] == 0

mij[entry_, matrix_] := 
 GetMatrixEntry[First[First[entry]], Last[First[entry]], matrix]

mkl[entry_, matrix_] := 
 GetMatrixEntry[First[Last[entry]], Last[Last[entry]], matrix]

muimuj[entry_, resultVect_] := 
 Extract[resultVect, First[First[entry]]]/
  Extract[resultVect, Last[First[entry]]]

mukmul[entry_, resultVect_] := 
 Extract[resultVect, First[Last[entry]]]/
  Extract[resultVect, Last[Last[entry]]]

auxUniqePairs[baseList_] := 
 Select[Tuples[baseList, 2], First[#] != Last[#] &]

COP2ViolationList[matrix_, resultList_] := 
 Select[Map[
   List[((mij[#, matrix] <= 
          mkl[#, matrix]) && (muimuj[#, resultList] <= 
          mukmul[#, resultList])) || ((mij[#, matrix] >= 
          mkl[#, matrix]) && (muimuj[#, resultList] >= 
          mukmul[#, resultList])), #] &, 
   auxUniqePairs[auxUniqePairs[Range[First[Dimensions[matrix]]]]]], 
  First[#] == False &]

COP2Check[matrix_, resultList_] := 
 Length[COP2ViolationList[matrix, resultList]] == 0

RankToVector[list_] :=
  Transpose[List[list]]

(***************** Discrepancy - Implementation *****************)

ErrorMatrix[M_, mju_] := 
  MapIndexed[1/#1 * Extract[mju, First[#2]]/Extract[mju, Last[#2]] &, M, {2}]

LocalDiscrepancyMatrix[M_, mju_] := 
  MapIndexed[Max[#1 - 1, 1/#1 - 1] & , ErrorMatrix[M, mju], {2}]

GlobalDiscrepancy[M_, mju_] := Max[Flatten[LocalDiscrepancyMatrix[M, mju]]]

RankOrder[rankList_] := Flatten[Map[Last[#] &,
	Reverse[SortBy[MapIndexed[{#1, #2} &, rankList],First[#1] &]]]]

End[]
EndPackage[]


PrincipalEigenValue[{{1, 2, 3}, {4, 5, 6}, {7, 3, 9}}]





mat = {{1, 2, 3}, {4, 5, 6}, {7, 3, 9}}
GeometricMean[mat]






