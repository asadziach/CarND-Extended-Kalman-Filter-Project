#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {

private:
	inline void update_(const Eigen::VectorXd &z,const Eigen::VectorXd &z_pred,
			Eigen::MatrixXd &H,Eigen::MatrixXd &R);
public:

	// state vector
	Eigen::VectorXd x_;

	// state covariance matrix
	Eigen::MatrixXd P_;
	/**
	 * Constructor
	 * @param x state vector
	 * @param P state ovariance matrix
	 */
	KalmanFilter();

	/**
	 * Destructor
	 */
	virtual ~KalmanFilter();

	/**
	 * Prediction Predicts the state and the state covariance
	 * using the process model
	 * @param delta_T Time between k and k+1 in s
	 * F = state transition matrix
	 * Q = process covariance matrix
	 */
	void Predict(Eigen::MatrixXd &F, Eigen::MatrixXd &Q);

	/**
	 * Updates the state by using standard Kalman Filter equations
	 * @param z The measurement at k+1
	 * @param H measurement matrix
	 * @param R measurement covariance matrix
	 */
	void Update(const Eigen::VectorXd &z, Eigen::MatrixXd &H,
			Eigen::MatrixXd &R);

	/**
	 * Updates the state by using Extended Kalman Filter equations
	 * @param z The measurement at k+1
	 * @param H measurement matrix
	 * @param R measurement covariance matrix
	 */
	void UpdateEKF(const Eigen::VectorXd &z, Eigen::MatrixXd &H,
			Eigen::MatrixXd &R);

	void DebugPrint();

};

#endif /* KALMAN_FILTER_H_ */
