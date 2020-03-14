#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

double const ESP = 1e-6;
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ += K*y;
  P_ = (I - K*H_)* P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  double h_rho = pow(px*px + py*py, 0.5);
  double h_phi = atan2(py, px);
  
  // check if h_rho is non-zero
  if (h_rho < ESP)
  {
    h_rho = ESP;
  }

  double h_rhodot = (px*vx + py*vy) / h_rho;

  VectorXd h = VectorXd(3);
  h << h_rho, h_phi, h_rhodot;

  VectorXd y = z - h;

  Tools tool = Tools();
  MatrixXd Hj = tool.CalculateJacobian(x_);

  MatrixXd S = Hj * P_ * Hj.transpose() + R_;
  MatrixXd K = P_* Hj * S.inverse();
  
  x_ += K * y;
  
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  P_  = (I - K * Hj) * P_;

}

