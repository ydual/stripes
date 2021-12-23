#include "Stripes.h"
#include <iostream>
#include <igl/readOBJ.h>

int main(void) {
    Eigen::MatrixXd V;
    Eigen::MatrixXd D;
    Eigen::MatrixXi F;
    Eigen::VectorXi branchIndex;
    Eigen::MatrixXd parameterization;
    Eigen::MatrixXi zeroIndex;
    std::vector<bool> isBorder;
    double theta = 0;
    double frequency = 130.;
    int numCoordinateFunctions = 1;
    igl::readOBJ("../data/square10.obj",V,F);


    D.resize(V.rows(),3);
    for(int i=0; i<V.rows(); ++i){
        D(i,0) = rand()/double(RAND_MAX);
        D(i,1) = rand()/double(RAND_MAX);
        D(i,2) = rand()/double(RAND_MAX);
    }


    int res = DDG::computeStripePatterns(V,F,D,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    
    return 0;
}
