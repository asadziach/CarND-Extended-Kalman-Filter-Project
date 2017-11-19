#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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

	/* The R matrix values and Q noise values are provided. Fine tuning
	 * is not required for this project. Unscented Kalman Filter project
	 * up next will explore this.
	 * */
	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0, 0, 0.0225;

	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

	/**_
	 * Finish initializing the FusionEKF.
	 * Set the process and measurement noises
	 */

	//the initial transition matrix F_
	F_ = MatrixXd(4, 4);
	F_ << 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1;

	//measurement matrix
	H_laser_ = MatrixXd(2, 4);
	H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;

	Q_ = MatrixXd(4, 4);
	Hj_ = MatrixXd(3, 4);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (!is_initialized_) {
		/**
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 */

		//state covariance matrix P
		ekf_.P_ = MatrixXd(4, 4);
		ekf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

		ekf_.x_ = VectorXd(4);

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			 */
			// Range - Radial distance from origin
			float rho = measurement_pack.raw_measurements_[0];
			// Bearing - angle between rho and x
			float phi = measurement_pack.raw_measurements_[1];
			// Radial Velocity - change of rho - range rate
			float rho_dot = measurement_pack.raw_measurements_[2];

			float px = rho * cos(phi);
			float py = rho * sin(phi);
			/* Although radar gives velocity data in the form of the range rate,
			 * a radar measurement does not contain enough information to determine
			 * the state variable velocities v​x and v​y. You can, however, use
			 * the radar measurements ρ and ϕ to initialize the state variable
			 * locations p​x and p​y.
			 * float vx = rho_dot * cos(phi);
			 * float vy = rho_dot * sin(phi);
			 * */
			ekf_.x_ << px, py, 0, 0;

		} else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			 */
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

		}

		previous_timestamp_ = measurement_pack.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	/**
	 * Update the state transition matrix F according to the new elapsed time.
	 - Time is measured in seconds.
	 * Update the process noise covariance matrix.
	 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	 */

	int noise_ax = 9;
	int noise_ay = 9;

	//compute the time elapsed between the current and previous measurements. dt - expressed in seconds
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	F_(0, 2) = dt;
	F_(1, 3) = dt;

	//set the process covariance matrix Q

	Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, 0, dt_4 / 4
			* noise_ay, 0, dt_3 / 2 * noise_ay, dt_3 / 2 * noise_ax, 0, dt_2
			* noise_ax, 0, 0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

	ekf_.Predict(F_, Q_);

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/
	/**
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	 */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.UpdateEKF(measurement_pack.raw_measurements_, Hj_, R_radar_);
	} else {
		// Laser updates
		ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);

	}

	ekf_.DebugPrint();
}
