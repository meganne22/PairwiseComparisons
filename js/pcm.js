/**
 * Created by Meg on 2015-01-12.
 */

var PCMethod = {};

/***********************************************************************/
/***********************Eigenvalue Rank Methods************************/
/**********************************************************************/

/*PrincipalEigenValue::usage = "PrincipalEigenValue[M] returns the largest principal eigenvalue of M";*/
PCMethod.PrincipalEigenValue = function(M) {
    return math.max(numeric.eig(M).lambda.x);
};

/*PrincipalEigenVector::usage = "PrincipalEigenVector[M] returns the eigenvector of M corresponding to its principal eigenvalue ";*/
PCMethod.PrincipalEigenVector = function(M) {
    var eigval = PCMethod.PrincipalEigenValue(M);
    var values = numeric.eig(M).E.x;
    var i =0;
    do {
        if (values[i] = eigval){
            return (numeric.eig(M).E.x[i]);
        }
        i++;
    } while (i <= values.length || values[i-1] == eigval);
};

/*SaatyIdx::usage = "SaatyIdx[M] returns the value of the Saaty Inconsistency Index computed for the matrix M";*/
PCMethod.SaatyIdx = function(M){
    var n = M.length;
    var alfa = PCMethod.PrincipalEigenValue(M);
    return (alfa -n)/(n -1);
};

/*EigenvalueRank::usage = "EigenvalueRank[M] returns the principal eigen vector of M rescalled in a way that the sum of its entries is 1";*/
PCMethod.EigenvalueRank = function(M){
    var A = PCMethod.PrincipalEigenVector(M);
    var sum = 0;
    for (var i=0; i< A.length; i++){
        sum += A[i];
    };
    for (var j=0; j< A.length; j++){
        A[j] = A[j]/sum;
    };
    return A;
};

/***********************************************************************/
/**************************Geometric Mean Rank Methods******************/
/***********************************************************************/

/*GeometricRank::usage = "GeometricRank[M] returns rank list given as geometric means of rows of the matrix M";*/
PCMethod.GeometricRank = function(M){
    var N = [];
    for (var i=0; i < M.length; i++){
        var sum = 1;

        for (var j=0; j < M.length; j++){
            sum= sum*M[i][j];
        };
        N[i] = Math.pow(sum,(1/M.length));
    };
    return N;
};

/*GeometricRescaledRank::usage = "GeometricRescaledRank[M] returns a geometric mean rank list of M rescalled in a way that the sum of its entries is 1";*/
PCMethod.GeometricRescaledRank = function(M){
    var A = PCMethod.GeometricRank(M);
    var sum = 0;
    for (var i=0; i< A.length; i++){
        sum += A[i];
    };
    for (var j=0; j< A.length; j++){
        A[j] = A[j]/sum;
    };
    return A;
};

/*******************************************************************************/
/****HRE (Heuristic Rating Estimation) - rankings with the reference values*****/
/******************************************************************************/
PCMethod.zeroIndices = function (v) {
    var z = [];
    var j = 0;
    for (var i = 0; i < v.length; i++){
        if (v[i] == 0){
            z.push(i);
        };
        j++;
    }
    return math.sort(z);
}
PCMethod.nonZeroIndices = function(v){
    var z = [];
    var j = 0;
    for (var i = 0; i < v.length; i++){
        if (v[i] != 0){
            z.push(i);
        };
        j++;
    }
    return math.sort(z);
}
/*HREConstantTermVector::usage = "HREConstantTermVectorEntry[M,v] returns the constant term vector for the HRE method. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREConstantTermVector = function(M,v){
    var zero = PCMethod.zeroIndices(v);
    var nonzero = PCMethod.nonZeroIndices(v);
    var termVec = [];
    var d = 1/ (M.length - 1);

    for (var n = 0; n < zero.length; n++){
        var elem = 0;
        for (var k = 0; k <nonzero.length; k++){
            elem += M[zero[n]][nonzero[k]] * v[k] ;
        };
        termVec[n] = elem * d;
    }
    return termVec;
};

/*HREMatrix::usage = "A = HREMatrix[M,v] returns a matrix for the HRE method. Togeter with b=HREConstantTermVector[M, v] it forms the linear equation system Au=b where u is the partial ranking vector. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREMatrix = function(M,v) {
    var A = [];
    while(A.push([]) < v.length-1);
    v = PCMethod.zeroIndices(v);

    for (var i=0; i < v.length; i++){
        for (var j = 0; j < v.length; j++){
            console.log(v);
            console.log(M[v[i]][v[j]]);
            A[i][j] = M[v[i]][v[j]];
        }
    };
    A = math.multiply(A, -1 / (M.length - 1));

    for (var n =0; n < A.length; n++){
        for (var m=0; m < A.length; m++){
            if (n == m){
                A[n][m] = 1;
            };
        };
    };

    return A;
};

/*HREPartialRank::usage = "HREPartialRank[M,v] returns values for the unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREPartialRank = function(M, v) {
    var termVector = PCMethod.HREConstantTermVector(M, v);
    var termM = math.ones(termVector.length)._data;

    for (var i=0; i <termVector.length; i++){
        termM[i] = termVector[i];
    };
    var HREMat = PCMethod.HREMatrix(M, v);
    var result = numeric.solve(HREMat, termM);

    return result;
};

/*HREFullRank::usage = "HREFullRank[M,v] returns values for both: known and unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREFullRank = function(M, v) {
    var partial = PCMethod.HREPartialRank(M, v);
    var i = 0;
    for (var j=0; j < v.length; j++){
        if (v[j] == 0){
            v[j] = partial[i];
            i++;
        }
    };
    return v;
};

/*HRERescaledRank::usage = "HRERescaledRank[M,v] returns the rescaled HREFullRank list, rescalled in a way that the entries sum up to 1"*/
PCMethod.HRERescaledRank = function(M,v) {
    var A = PCMethod.HREFullRank(M,v);
    var sum = 0;
    for (var i=0; i< A.length; i++){
        sum += A[i];
    };
    for (var j=0; j< A.length; j++){
        A[j] = A[j]/sum;
    };
    return A;
};

/*******************************************************************************/
/****HRE Geometric (Heuristic Rating Estimation using Geometric Mean) - ********/
/***********************rankings with the reference values**********************/
/******************************************************************************/

/*HREGeomConstantTermVector::usage = "HREGeomConstantTermVectorEntry[M,v] returns the constant term vector for the HRE method. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREGeomConstantTermVector = function(M,v){
    var zero = PCMethod.zeroIndices(v);
    var nonzero = PCMethod.nonZeroIndices(v);
    var termVec = [];
    var d = 1/ (M.length - 1);

    if (nonzero.length ==0){
        return termVec;
    }

    var multiplier =  1;
    for(var i=0; i<nonzero.length; i++){
        multiplier *= v[nonzero[i]];
    }

    for (var n = 0; n < zero.length; n++){
        var elem = 1;
        for (var k = 0; k <M.length; k++){
            elem *= M[zero[n]][k] ;
        };
        termVec.push(math.log10(elem * multiplier));
    }
    return termVec;
};

 /*HREGeomMatrix::usage = "A = HREGeomMatrix[M,v] returns a matrix for the HRE Geometric method. Togeter with b=HREGeomConstantTermVector[M, v] it forms the linear equation system Au=b where u is the intermediate partial ranking vector. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREGeomMatrix = function(M,v){
    var zeros = PCMethod.zeroIndices(v);
    var result = [];
    while(result.push([]) < zeros.length);

    for (var i=0; i<result.length; i++){
        for(var j=0; j<result.length; j++){
            if (i ==j){
                result[i][j] = M.length - 1;
            }
            else{
                result[i][j] = -1;
            }

        }
    };
    return result;
};

 /*HREGeomIntermediateRank::usage = "HREGeomIntermediateRank[M,v] returns values for the unknown alternatives before they are raised to the power 10. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREGeomIntermediateRank = function(M,v){
    var termVector = PCMethod.HREGeomConstantTermVector(M,v);
    var termVec = [] ;

    for(var i=0; i<termVector.length; i++){
        termVec[i] = termVector[i];
    }
    var HREMat = PCMethod.HREGeomMatrix(M,v);
    var result = numeric.solve(HREMat, termVec);

    return result;
};

 /*HREGeomPartialRank::usage = "HREGeomPartialRank[M,v] returns values for the unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREGeomPartialRank = function(M,v){
    var intermediateRank = PCMethod.HREGeomIntermediateRank(M,v);
    var RankLinked = intermediateRank;
    var tpReplace = 0;

    for(var i=0; i<RankLinked.length; i++){
        toReplace = math.pow(10,RankLinked[i]);
        RankLinked[i] = toReplace;
    };
    return RankLinked;
};

 /*HREGeomFullRank::usage = "HREGeomFullRank[M,v] returns values for both: known and unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREGeomFullRank = function(M,v){
    var partial = PCMethod.HREGeomPartialRank(M,v);
    var i =0;
    for (var j=0; j< v.length; j++){
        if (v[j] == 0) {
            v[j] = partial[i];
            i++;
        }
    }
    return v;
};

 /*HREGeomRescaledRank::usage = "HREGeomRescaledRank[M,v] returns the rescaled HREGeomFullRank list, rescalled in a way that the entries sum up to 1"*/
PCMethod.HREGeomRescaledRank = function(M,v){
    var vector = PCMethod.HREGeomFullRank(M,v);
    var sum = 0;
    for (var i=0; i<vector.length; i++){
        sum +=vector[i];
    }
    for (var j=0; j<vector.length; j++){
        vector[j] = vector[j]/sum;
    }
    return vector;
};

/*******************************************************************************/
/***********************Koczkodaj Inconsistency Methods**********************/
/******************************************************************************/

/*KoczkodajTriadIdx::usage = "KoczkodajTriadIdx[{a,b,c}] returns the Koczkodaj triad inconsistency i.e. Min[Abs[1-c/ab],Abs[1-ab/c]]";
KoczkodajIdx::usage = "KoczkodajIdx[M] returns the value of the Koczkodaj inconsistency computed for the matrix M";
KoczkodajTheWorstTriad::usage = "KoczkodajTheWorstTriad[M] returns the worst triad in M according to Koczkodaj inconsistency criterion";
KoczkodajTheWorstTriads::usage = "KoczkodajTheWorstTriads[M, n] returns n worst triads in M according to Koczkodaj inconsistency criterion";
KoczkodajConsistentTriad::usage = "KoczkodajConsistentTriad[{a,b,c}] returns the closest consistent triad";
KoczkodajImproveMatrixStep::usage = "KoczkodajImproveMatrixStep[M] returns improved (more consistent) matrix in which the most inconsistent triad is replaced by the closest consistent one";
*/


/*******************************************************************************/
/***********Condition Of Order Preservation by Bana e Costa and Vansnick************/
/******************************************************************************/

/*COP1ViolationList::usage = "COP1ViolationList[M,r] returns the list of indices that violate the first Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick). Here, M is the PC Matrix, whilst r is the ranking result list."

 COP1Check::usage = "COP2ViolationList[M,r] it is true when the first Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick) is violated. Here, M is the PC Matrix, whilst r is the ranking result list."

 COP2ViolationList::usage = "COP2ViolationList[M,r] returns the list of indices that violate the second Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick). Here, M is the PC Matrix, whilst r is the ranking result list."

 COP2Check::usage = "COP2ViolationList[M,r] it is true when the second Condition of Order Preservation (COP) postulate (formulated by Bana e Costa and Vansnick) is violated. Here, M is the PC Matrix, whilst r is the ranking result list."*/

/*zwraca liste indeksow gdzie wystepuja rozbieznosci*/
PCMethod.COP1Violation = function(M, r){
 var violations = [];
    for(var i=0; i< M.length; i++){
        for(var j=0; j< M.length; j++){
            if(i != j){
                if(M[i][j] > 1){
                    if(!(r[i] > r[j])){
                        violations.push([i,j]);
                    }
                }
            }
        }
    }
    return violations;
};

PCMethod.COP2Violation= function(M, r){
    var violations = [];

    for(var i=0; i< M.length; i++){
        for(var j=0; j< M.length; j++){
            for(var k=0; k< M.length; k++){
                for(var l=0; l< M.length; l++){
                    if(i != j && k != l){
                        if(M[i][j] > M[k][l]) {
                            if (!(r[i]/ r[j] > r[k]/r[l])) {
                                violations.push([i, j, k, l]);
                            }
                        }
                    }
                }
            }
        }
    }
    return violations;
};
/*******************************************************************************/
/**************************Ranking Discrepancy*******************************/
/******************************************************************************/


/*ErrorMatrix::usage = "ErrorMatrix[M, mju] computes matrix containing entries in form e_ij = m_ji*mju(c_i)/mju(c_j). When the matrix is consistent every e_ij equals 1"

 LocalDiscrepancyMatrix::usage = "LocalDiscrepancyMatrix[M, mju] computes matrix with entries d_ij = max{e_ij - 1, 1/e_ij - 1}"

 GlobalDiscrepancy::usage = "GlobalDiscrepancy[M, mju] returns maximal entry of  LocalDiscrepancyMatrix[M, mju]"*/

PCMethod.ErrorMatrix = function(){};



/*******************************************************************************/
/*************************Auxiliary functions*******************************/
/******************************************************************************/

/*GetMatrixEntry[i ,j, M] - the same as M[i,j]*/
PCMethod.GetMatrixEntry = function(x, y, M){
    var v = M[x][y];
    console.log(v);
    return v;
};

/*RecreatePCMatrix::usage = "RecreatePCMatrix[M] recreates reciprocal matrix on the basis of upper-triagle of M";*/
PCMethod.RecreatePCMatrix = function (M) {
    for (var i = 0; i < M.length; i++){
        for(var j = 0 ; j < M.length; j++){
            if (i < j){
                M[i][j] = M[i][j];
            }
            else if (i == j){
                M[i][j] = 1;
            }
            else {
                M[i][j] = math.divide(1, PCMethod.GetMatrixEntry(j,i,M));
            }
        }
    }
    return M;
};