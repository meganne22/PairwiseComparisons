/**
 * Created by Meg on 2015-01-12.
 */
var math = require('E:\\Git\\PairwiseComparisons\\Projekt\\mathjs-master\\index');
//Eigenvalue Rank Methods
//
//PrincipalEigenValue::usage = "PrincipalEigenValue[M] returns the largest principal eigenvalue of M";
//PrincipalEigenVector::usage = "PrincipalEigenVector[M] returns the eigenvector of M corresponding to its principal eigenvalue ";
//SaatyIdx::usage = "SaatyIdx[M] returns the value of the Saaty Inconsistency Index computed for the matrix M";
//EigenvalueRank::usage = "EigenvalueRank[M] returns the principal eigen vector of M rescalled in a way that the sum of its entries is 1";

var PCMethod = {};

PCMethod.PrincipalEigenValue = function(M) {
    return max(math.eval(M));
};