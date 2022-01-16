#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<math.h>
#include"functions.h"

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd function_y(MatrixXd para, MatrixXd x) {
	int number_of_data_points = x.rows();
	MatrixXd y(number_of_data_points, 1);
	for (int i = 0; i < number_of_data_points; i++) {

		//y(i, 0) = para(0, 0) * pow(x(i, 0), 2) + para(1, 0) * exp(pow(z(i, 0), 2) + 1) + para(2, 0);

		y(i, 0) = para(0, 0) * pow(x(i, 0), 2) + para(1, 0) * exp(para(2, 0) * x(i, 0)); // +para(3, 0) * x(i, 0);
	}
	return y;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

MatrixXd getJacobianMatrix(MatrixXd para_est, MatrixXd deflection, MatrixXd ym, MatrixXd input) {

	MatrixXd Jacobian_Matrix(ym.rows(), para_est.rows());
	MatrixXd y = function_y(para_est, input);
	MatrixXd y_deflected(ym.rows(), 1);

	for (int i = 0; i < para_est.rows(); i++) {

		para_est(i, 0) = para_est(i, 0) + deflection(i, 0);		/*Changing the parameters one by one */

		y_deflected = function_y(para_est, input);				/*Computing the deflected function arrray */
		for (int j = 0; j < input.rows(); j++) {

			// [f(v, p + dp) - f(v, p) ] / [dp] 

			Jacobian_Matrix(j, i) = (y_deflected(j, 0) - y(j, 0)) / deflection(i, 0);
		}
		para_est(i, 0) = para_est(i, 0) - deflection(i, 0);		/*Bringing back the parametes to original value*/
	}
	return Jacobian_Matrix;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/


MatrixXd GradientDescent(MatrixXd para_guess, MatrixXd deflection, MatrixXd ym, MatrixXd input) {

	cout << "-> Entered Gradient Descent\n";

	//cout << "para_guess: " << para_guess << endl;
	/*cout << "deflection: " << deflection << endl;
	cout << "ym: " << ym << endl;
	cout << "input: " << input << endl;*/

	int npara = para_guess.rows(), ndata = ym.rows();

	//MatrixXd IdentityMat = MatrixXd::Identity(npara, npara);
	//MatrixXd H(npara, npara);
	MatrixXd d(ndata, 1);
	MatrixXd J(ndata, npara);
	double error,error_gd;

	MatrixXd y_init = function_y(para_guess, input);
	//cout << "yinit: " << y_init << endl;

	//double lambda = 10;
	//int updateJ = 1;
	double alpha = 1e-6;
	MatrixXd para_est = para_guess;
	MatrixXd para_gd, y_est_gd, y_est;
	int maxiter = 10000, counter=0;
	

	while (counter < maxiter) {
		
		cout << "--> Iteration : " << counter << endl;

		if (counter == 0) {
			y_est = function_y(para_est, input);
			d = ym - y_est;
		}
		else {
			y_est = y_est_gd;
			para_est = para_gd;
			error = error_gd; 
		}

		J = getJacobianMatrix(para_est, deflection, ym,input);
		//cout << J << endl;
		MatrixXd dp = alpha * J.transpose() * d;
		cout << "dp: \n" << dp;
		para_gd = para_est + dp;
		y_est_gd= function_y(para_gd, input);
		MatrixXd d_gd = ym - y_est_gd;
		MatrixXd temp2 = d_gd.transpose() * d_gd;
		error_gd = temp2(0,0);
		
		if (dp.norm() < 1e-8) {
			counter = maxiter;
		}
		else counter++;
		
		
	}

	return para_gd;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
