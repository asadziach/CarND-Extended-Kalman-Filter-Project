#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict(Eigen::MatrixXd &F, Eigen::MatrixXd &Q) {
  x_ = F * x_;
  MatrixXd Ft = F.transpose();
  P_ = F * P_ * Ft + Q;
}

// update the state by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {

  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
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

void KalmanFilter::UpdateEKF(const VectorXd &z, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}

void KalmanFilter::DebugPrint(){

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}
