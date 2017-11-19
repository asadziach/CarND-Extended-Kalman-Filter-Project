#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {
}

void KalmanFilter::Predict(Eigen::MatrixXd &F, Eigen::MatrixXd &Q) {
	x_ = F * x_;
	MatrixXd Ft = F.transpose();
	P_ = F * P_ * Ft + Q;
}

void KalmanFilter::update_(const Eigen::VectorXd &z,const Eigen::VectorXd &z_pred,
		Eigen::MatrixXd &H,Eigen::MatrixXd &R){

	VectorXd y = z - z_pred;

	// Normalize angle if result is not in the range pi to -pi
	y[1] = atan2(sin(y[1]), cos(y[1]));

	MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H) * P_;
}

// update the state by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z, Eigen::MatrixXd &H,
		Eigen::MatrixXd &R) {

	VectorXd z_pred = H * x_;

	update_(z,z_pred,H,R);
}

/**
 * update the state by using Extended Kalman Filter equations
 */
void KalmanFilter::UpdateEKF(const VectorXd &z, Eigen::MatrixXd &H,
		Eigen::MatrixXd &R) {
	/*
	 * convert x to polar coordinates, the function h(x) maps values from Cartesian coordinates
	 * to polar coordinates
	 */
	float rho = sqrt((x_(0) * x_(0)) + (x_(1) * x_(1)));

	/* atan2 internally checks for divide by zero and returns angle between pi and -pi.
	 * no need to normalize.
	 */
	float phi = atan2(x_(1), x_(0));

	float rhodot = 0;
	//check division by zero
	if (fabs(rho) > 0.0001) {
		rhodot = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
	}

	VectorXd z_pred(3);
	z_pred << rho, phi, rhodot;

	update_(z,z_pred,H,R);
}

void KalmanFilter::DebugPrint() {

	// print the output
	cout << "x_ = " << x_ << endl;
	cout << "P_ = " << P_ << endl;

}
