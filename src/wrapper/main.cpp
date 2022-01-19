#include "Stripes.h"
#include <iostream>
#include <igl/readOBJ.h>

double computeObjective(Eigen::VectorXd& param_cur, Eigen::VectorXd& param_tar){
	double obj = 0;
	for(int i=0; i<param_cur.rows(); i++){
		obj+=pow((param_cur(i)-param_tar(i)),2);
	}
	return obj;
}


int main(void) {
    Eigen::MatrixXd V;
    Eigen::MatrixXd D_cur, D_tar;
    Eigen::MatrixXi F;
    Eigen::VectorXi branchIndex;
    Eigen::MatrixXd parameterization;
    Eigen::MatrixXi zeroIndex;
    std::vector<bool> isBorder;
    double theta = 0;
    double frequency = 130.;
    int numCoordinateFunctions = 1;
    igl::readOBJ("../data/square10.obj",V,F);

    srand(20);

    D_cur.resize(V.rows(),3);
    // for(int i=0; i<V.rows(); ++i){
    //     D_cur(i,0) = 0.6+0.1*(i*i*i);
    //     D_cur(i,1) = 0.8-0.8*(i*i);
    //     D_cur(i,2) = 0;
    //     D_cur.row(i).normalize();
    // }
    for(int i=0; i<V.rows(); ++i){
        D_cur(i,0) = (rand()+1.)/(RAND_MAX+1.);
        D_cur(i,1) = (rand()+1.)/(RAND_MAX+1.);
        D_cur(i,2) = 0;
        //D_cur.row(i).normalize();
    }

    std::cout<<D_cur.transpose()<<std::endl;

    // D_tar.resize(V.rows(),3);
    // for(int i=0; i<V.rows(); ++i){
    //     D_tar(i,0) = D_cur(i,0)-0.03;
    //     D_tar(i,1) = D_cur(i,1)+0.03;
    //     D_tar(i,2) = D_cur(i,1)-0.03;
    // }

    Eigen::VectorXd param_cur, param_tar, param_p, param_pp, v_tar;
    std::vector<Eigen::MatrixXd> dpdD, dummy;

    Eigen::MatrixXd D_init = D_cur;

    param_cur.setZero();
    param_tar.setZero();
    param_p.setZero();
    param_pp.setZero();
    v_tar.setZero();

    //for(int i=0; i<3; ++i)
        DDG::computeStripePatterns(1,v_tar,dpdD,param_cur,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    //DDG::computeStripePatterns(0,v_tar,dpdD,param_cur,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    //std::cout<<dpdD[0].rows()<<std::endl;

    // double eps = 1e-5;
    // for(int i=0; i<D_cur.rows(); ++i)
    // {
    //     D_cur(i,1) += eps;
    //     for(int j=0; j<3; ++j)
    //         DDG::computeStripePatterns(0,v_tar,dummy,param_p,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    //     D_cur(i,1) -= 2*eps;
    //     for(int j=0; j<3; ++j)
    //         DDG::computeStripePatterns(0,v_tar,dummy,param_pp,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    //     D_cur(i,1) += eps;

    //     // std::cout<<param_cur.transpose()<<std::endl;
    //     // std::cout<<param_p.transpose()<<std::endl;

    //     for(int j=0; j<param_cur.rows(); ++j){
    //         double dpdDN = (param_p(j)-param_pp(j))/(2*eps);
    //         if(abs(dpdDN) > 1e-7)
    //             std::cout<<i<<" "<<j<<" "<<dpdDN<<" "<<dpdD[i](j,1)<<" "<<(dpdDN-dpdD[i](j,1))/dpdDN<<std::endl;
    //     }   
        
    // }

    // DDG::computeStripePatterns(1,v_tar,dpdD,param_cur,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    // DDG::computeStripePatterns(0,v_tar,dummy,param_tar,V,F,D_tar,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    // std::cout<<computeObjective(param_cur, param_tar)<<std::endl;`   
    // // std::cout<<D_cur.transpose()<<std::endl;

    // int iter = 0;
    // double prev = computeObjective(param_cur, param_tar);

    // while(iter<1000)
    // {
    //     double alpha_init = 0.00001;
	// 	double tao = 0.5;
	// 	double beta = 0.001;
    //     double prev_obj = computeObjective(param_cur, param_tar);


    //     while(true)
    //     {
    //         Eigen::MatrixXd dOdDs(V.rows(),3);
    //         dOdDs.setZero();

    //         double gdp = 0;

    //         for(int i=0; i<V.rows(); ++i){
    //             for(int j=0; j<3; ++j){
    //                 for(int k=0; k<2*V.rows(); ++k){
    //                     if(dpdD[i](k,j) > 1e-10)
    //                         dOdDs(i,j) += (param_cur(k)-param_tar(k))*dpdD[i](k,j);
    //                     gdp -= dOdDs(i,j)*dOdDs(i,j);
    //                     D_cur(i,j) -= alpha_init * dOdDs(i,j);
    //                 }
    //             }
    //         }
    //         //std::cout<<D_cur.transpose()<<std::endl;
    //         DDG::computeStripePatterns(0,v_tar,dpdD,param_cur,V,F,D_cur,theta,frequency,numCoordinateFunctions,branchIndex,parameterization,zeroIndex,isBorder);
    //         alpha_init*=tao;
    //         // std::cout<<"--------------------------------------------------"<<std::endl;
	// 		// std::cout<<"Two values: "<<computeObjective(param_cur, param_tar)<<" "<<prev_obj+alpha_init*beta*gdp<<std::endl;
	// 		// std::cout<<"--------------------------------------------------"<<std::endl;
    //         if(computeObjective(param_cur, param_tar)<= prev_obj+alpha_init*beta*gdp) break;
    //         prev_obj = computeObjective(param_cur, param_tar);
    //     }

    //     // if(prev < computeObjective(param_cur, param_tar)) break;
    //     // prev = computeObjective(param_cur, param_tar);

    //     //std::cout<<D_cur.transpose()<<std::endl;
    //     std::cout<<computeObjective(param_cur, param_tar)<<std::endl;
    //     iter++;
    // }

    // std::cout<<D_cur<<std::endl;

    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file("../data/square200.csv");
    if (file.is_open())
    {
        file << D_cur.format(CSVFormat);
        file.close();
    }

    // std::ofstream file1("../data/square10_tar.csv");
    // if (file1.is_open())
    // {
    //     file1 << D_tar.format(CSVFormat);
    //     file1.close();
    // }

    // std::ofstream file2("../data/square10_cur.csv");
    // if (file2.is_open())
    // {
    //     file2 << D_cur.format(CSVFormat);
    //     file2.close();
    // }
    
    return 0;
}
