#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  P_init_ = MatrixXd(4,4);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0; 

  P_init_ << 10,  0,  0,  0,
             0, 10,  0,  0,
             0,  0, 1000, 0,
             0, 0, 0, 1000;

  //Hj_ = Tools().CalculateJacobian();
  ekf_ = KalmanFilter();

  //have not set the process and measurement noises 
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

MatrixXd FusionEKF::CalculateProcessCovarianceMatrix(double dt, double sigma_x_square, double sigma_y_square)
{
  MatrixXd Q = MatrixXd(4,4);
  double dt_2 = pow(dt, 2);
  double dt_3 = pow(dt, 3);
  double dt_4 = pow(dt, 4);
  
  Q << dt_4 * sigma_x_square / 4,                  0,  dt_3 * sigma_x_square / 2,       0,
                              0, dt_4 * sigma_y_square / 4,                   0, dt_3 * sigma_y_square / 2,
                              dt_3 * sigma_x_square / 2, 0,       dt_2 * sigma_x_square,      0,
                              0, dt_3 * sigma_y_square / 2,                   0, dt_2 * sigma_y_square;        

  return Q;
}

MatrixXd FusionEKF::CalculateTransitionMatrix(double dt)
{
  /*
  * Compute F - transition matrix 
  * @params: 
  */
  MatrixXd F = MatrixXd(4,4);
  F << 1, 0, dt, 0,
          0, 1, 0, dt,
          0, 0, 1, 0,
          0, 0, 0, 1;
  return F;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  double noise_ax = 9.0; // represent sigma x square
  double noise_ay = 9.0; // represent sigma y square
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    MatrixXd F_init = CalculateTransitionMatrix(0.0);
    MatrixXd Q_init = CalculateProcessCovarianceMatrix(0.0, noise_ax, noise_ay);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
    double rho = measurement_pack.raw_measurements_(0);
    double theta = measurement_pack.raw_measurements_(1);
    double rho_dot = measurement_pack.raw_measurements_(2);

    double px = rho * cos(theta);
    double py = rho * sin(theta);
    double vx = rho_dot * cos(theta);
    double vy = rho_dot * sin(theta);

    VectorXd x_in = VectorXd(4);
    x_in << px, py, vx, vy;
    Hj_ = tools.CalculateJacobian(x_in);

    ekf_.Init(x_in, P_init_, F_init, Hj_, R_radar_, Q_init);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      double px = measurement_pack.raw_measurements_(0);
      double py = measurement_pack.raw_measurements_(1);
      VectorXd x_in = VectorXd(4);
      x_in << px, py, 0, 0;
      
      ekf_.Init(x_in, P_init_, F_init, H_laser_, R_laser_, Q_init);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.Q_ = CalculateProcessCovarianceMatrix(dt, noise_ax, noise_ay);
  ekf_.F_ = CalculateTransitionMatrix(dt);

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
