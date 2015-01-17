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
PCMethod.EigenvalRank = function(M){
    var eigvec = PCMethod.PrincipalEigenVector(M);
    var sum = math.sum(eigvec);
    for (var i=0; i < eigvec.length; i++ ){
        eigvec[i] = eigvec[i]/sum;
    }
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
            console.log("wew petla: "+ sum);
        };
        N[i] = Math.pow(sum,(1/M.length));
        console.log("pierwsza petla " + sum + ", i: " + i + ", N: " +N);
    };
    return N;
};

/*GeometricRescaledRank::usage = "GeometricRescaledRank[M] returns a geometric mean rank list of M rescalled in a way that the sum of its entries is 1";*/
PCMethod.GeometricRescaledRank = function(){

};

/*******************************************************************************/
/****HRE (Heuristic Rating Estimation) - rankings with the reference values*****/
/******************************************************************************/

/*HREConstantTermVector::usage = "HREConstantTermVectorEntry[M,v] returns the constant term vector for the HRE method. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREConstantTermVector = function(){

};

/*HREMatrix::usage = "A = HREMatrix[M,v] returns a matrix for the HRE method. Togeter with b=HREConstantTermVector[M, v] it forms the linear equation system Au=b where u is the partial ranking vector. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREMatrix = function() {

};

/*HREPartialRank::usage = "HREPartialRank[M,v] returns values for the unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREPartialRank = function() {

};

/*HREFullRank::usage = "HREFullRank[M,v] returns values for both: known and unknown alternatives. M is the PC Matrix, whilst v is the vector of the known alternatives. In v the unknown alternatives are marked by the value 0"*/
PCMethod.HREFullRank = function() {

};

/*HRERescaledRank::usage = "HRERescaledRank[M,v] returns the rescaled HREFullRank list, rescalled in a way that the entries sum up to 1"*/
PCMethod.HRERescaledRank = function() {

};






/*RankOrder::usage = "RankOrder[rankList] - prints order of elements on the rank list"*/
PCMethod.RankOrder = function(){

};
