#include "kalman_filter.h"
#include <iostream>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  //KF Variables
  x_ = x_in; //object state
  P_ = P_in; //object covariance matrix
  F_ = F_in; //state transition matrix
  H_ = H_in; //measurement matrix
  R_ = R_in; //measurement covariance
  Q_ = Q_in; // process covariance
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
/*
     * KF Measurement update step
     */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  cout << "UpdateEKF.........." << endl;


  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float phi = 0;
  float rho_dot = 0;

  float rho = sqrt( px * px + py * py);

  // avoid division by zero
  if(fabs(px) < 0.0001){
    cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
  }
  else {
    phi = atan2(py, px);
  }

  // avoid division by zero
  if (rho < 0.0001) {
    rho = 0.0001;
    cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
  }
  else {
    rho_dot = (px*vx + py*vy) / rho;
  }

  VectorXd h = VectorXd(3);
  h << rho, phi, rho_dot; // For radar H * x becomes h(x)

  VectorXd y = z - h; // Using h instead of Jacobian Hj_ here!



  while (y(1)>M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1)<-M_PI) {
    y(1) += 2 * M_PI;
  }


  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

