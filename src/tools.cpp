#include "tools.h"
#include <iostream>
#include <cmath>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout; 

const double ESP = 1e-6;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
     std::cout << "Invalid estimation or ground-truth data.\n";
     return rmse;
  }
   
   VectorXd accum(4);
   accum << 0, 0, 0, 0; 
   size_t n = estimations.size();
   
   for (size_t i = 0; i < n; ++i)
   {
      VectorXd residual = estimations[i] - ground_truth[i];
      VectorXd residual_square = residual.array() * residual.array();

      accum += residual_square;
   }

   rmse = accum / n;
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  double distance_square = px*px + py*py;

  
  // check if distance is zero.
  if (distance_square < ESP)
  {
     cout << "Division by zero. Setting distance to be " << ESP << ".\n";
     distance_square = ESP;

  }
   double distance = pow(distance_square, 0.5);
   Hj << px / distance,           py / distance,             0,             0,
 -py / distance_square,     px / distance_square,      0,             0,
   py * (vx*py - vy*px) / pow(distance, 3), px * (vy*px - vx*py)/ pow(distance, 3),    px / distance, py / distance;

  return Hj;
}
