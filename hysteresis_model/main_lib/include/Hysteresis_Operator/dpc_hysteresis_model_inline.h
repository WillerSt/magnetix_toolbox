/*
 * dpc_hysteresis_model_inline.h
 *
 *  Created on: 30.03.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_INLINE_H_
#define HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_INLINE_H_

template <int dimension, typename CS>
inline Eigen::Matrix<double, dimension, 1> DPC_hysteresis_model<dimension, CS>::get_H() const
{
	return this->Hnew;
}

template <int dimension, typename CS>
inline Eigen::Matrix<double, dimension, 1> DPC_hysteresis_model<dimension, CS>::get_M() const
{
	return this->value * hystLib::constants::mu0inv;
}

template <int dimension, typename CS>
inline Eigen::Matrix<double, dimension, 1> DPC_hysteresis_model<dimension, CS>::get_B() const
{
	return this->value + hystLib::constants::mu0 * this->Hnew;
}

template <int dimension, typename CS>
inline Eigen::Matrix<double, dimension, 1> DPC_hysteresis_model<dimension, CS>::get_J() const
{
	return this->value;
}

template <int dimension, typename CS>
inline int DPC_hysteresis_model<dimension, CS>::get_number_of_hysterons()
{
	return this->hysterons.size();
}


#endif /* HYSTERESIS_OPERATOR_DPC_HYSTERESIS_MODEL_INLINE_H_ */
