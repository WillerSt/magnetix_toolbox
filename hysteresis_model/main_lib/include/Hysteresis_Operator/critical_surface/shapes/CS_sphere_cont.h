/*
 * CS_sphere_cont.h
 *
 *  Created on: 27.11.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_CONT_H_
#define HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_CONT_H_

#include "Hysteresis_Operator/critical_surface/shapes/CS_sphere.h"
#include "hystLib_constants.h"
#include "hystLib_functions.h"
#include "Hysteresis_Operator/Misc_Tools/hyst_config.h"
#include <vector>
#include "Hysteresis_Operator/Misc_Tools/cont_int_table.h"
#include <Eigen/StdVector>

//#define CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS

//#define CS_SPHERE_CONT_USE_IDEALIZED_INITIALIZATION

#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
#include <iostream>
#include <iomanip>
#endif

#define RAD 50


//#define DPC_CONT_USE_LOOKUP // you shouldn't

template <int dimension>
class CS_sphere_cont : public CS_sphere<dimension>
{
private:
	const double DISTANCE_EPS = 1e-14;
	const double EQUALITY_EPS = 1e-14;
	const double INIT_DISTANCE_EPS = 1e-10;
public:
	using EigenType =  Eigen::Matrix<double, dimension, 1>;
	using EigenJacob = Eigen::Matrix<double, dimension, dimension>;
	using EigenStdVector = std::vector<EigenType,Eigen::aligned_allocator<EigenType>>;

public:
	const double radInner, radOuter, diffRad, diffRadInv, radInnerSq, radOuterSq;
	static const char shapeIdentifier = 'S';

public:
	// constructors
	CS_sphere_cont (const EigenType & positionIn, const double& radIn, const double& densityIn, const unsigned int& surfaceIDIn)
	:CS_sphere<dimension>(positionIn, radIn, densityIn, surfaceIDIn),
	 radInner(get_min_rad(radIn)), radOuter(get_max_rad(radIn)), diffRad(radOuter-radInner), diffRadInv(1/diffRad),
	 radInnerSq(radInner*radInner), radOuterSq(radOuter*radOuter)
	 {

	 }

	CS_sphere_cont (const EigenType & positionIn, const double& radIn, double radInnerIn, double radOuterIn,const double& densityIn, const unsigned int& surfaceIDIn)
	:CS_sphere<dimension>(positionIn, radIn,  densityIn, surfaceIDIn),
	 radInner(radInnerIn), radOuter(radOuterIn), diffRad(radOuter-radInner), diffRadInv(1/diffRad),
	 radInnerSq(radInner*radInner), radOuterSq(radOuter*radOuter)
	 {

	 }

private:
	double get_min_rad(const double& radIn) const
	{
		if ((radIn-RAD)<0.5||true)
		{
			return 0.8 * radIn;
		}
		else
		{
			return radIn -RAD;
		}
	}
	double get_max_rad(const double& radIn) const
	{
		if ((radIn-RAD)<0.5||true)
		{
			return 1.2 * radIn;
		}
		else
		{
			return radIn  + RAD;
		}
	}

public:

	double get_max_radius() const
	{
		return this->radOuter;
	}

	// calculate magnetization contributed by the whole path
	void calculate_frozen_part(const EigenStdVector & relevantPoints, const EigenType& direction,
				const std::vector<double>& distances, const EigenType& Hnew, EigenType& frozenPart, EigenType& unfrozenPart) const
		{


			// determine innermost and outermost relevant point
			const int sizeOfVectors = distances.size()-1;
			int indexOuter = sizeOfVectors;
			int indexInner = sizeOfVectors;
			for (int i = sizeOfVectors; i >= 0; --i)
			{
				if (distances[i] < radInner)
				{
					indexInner = i;
				}

				if (distances[i] > radOuter)
				{
					indexOuter = i;
					break;
				}
				else
				{
					indexOuter = i;
				}
			}


			#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
			if (indexOuter > indexInner)
			{
				std::cout << "WARNING: Are you sure about the indices?\n";
			}
			#endif

			//unfrozenPart = EigenType::Zero();

			EigenType frozenTemp = EigenType::Zero();

			// if the hysteron has never been fully unfrozen, assume demagnetized state (sum over all hysterons)
			if (indexOuter == indexInner)
			{
				if ((distances[indexInner] >radInner))
				{
					unfrozenPart += (distances[indexInner] - radInner) * diffRadInv * this->density * direction;
					frozenTemp += calculate_zero_state(distances[indexInner]) * this->density;
				}
				else
				{
					frozenTemp += calculate_zero_state(radInner)* this->density;
				}
				frozenPart += frozenTemp;

				#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
				if (std::isnan(frozenPart[0]) || std::isnan(frozenPart[1]) )
				{
					std::cout << "ERROR: NaN value in partial calculation of frozen value, initialization case" << std::endl;
					exit(1);
				}
				#endif
				return;
			}

			int calcCase = 4;

			for (int i = indexInner; i > indexOuter; i--)
			{
				calcCase = 0;
				if ( (i == indexInner) && ((i-1) != indexOuter))
				{
					if (distances[i] < radInner) // evolution ends outside of the ring
					{
						calcCase = 1;
					}
					else	// segment within ring
					{
						calcCase = 3;
						unfrozenPart += (distances[indexInner] - radInner) * diffRadInv * direction * this->density;
					}
				}
				else if((i == indexInner) && ((i-1) == indexOuter))
				{
					if ((distances[indexOuter] >= radOuter) && (distances[indexInner] < radInner)) // evolution passing surface in one step
					{
						calcCase = 5;

					}
					else if ((distances[indexOuter] >= radOuter) && (distances[indexInner] > radInner)) // partly unfrozen
					{
						calcCase = 2;
						unfrozenPart += (distances[indexInner] - radInner) * diffRadInv * direction * this->density ;
					}
					else if ((distances[indexOuter] < radOuter) && (distances[indexInner] < radInner)) // evolution started in ring
					{
						frozenTemp += calculate_zero_state(distances[indexOuter]);
						//frozenTemp +=  -(radOuter-distances[indexOuter]) * diffRadInv * this->position.normalized();
						calcCase = 1;
					}
					else
					{
						calcCase = 3;
						unfrozenPart += (distances[indexInner] - radInner) * diffRadInv * direction *this->density;
						frozenTemp +=  calculate_zero_state(distances[indexOuter]);
					}
				}
				else if ((i != indexInner) && ((i-1) == indexOuter))
				{
					if (distances[indexOuter] > radOuter) // evolution passing outer surface
					{
						calcCase = 2;
					}
					else //partly unfrozen
					{
						frozenTemp  += calculate_zero_state(distances[indexOuter]);
						calcCase = 3;
					}
				}
				else // segment is within the ring
				{
					calcCase = 3;
				}

				if (abs(distances[i-1] -distances[i]) < EQUALITY_EPS)
				{
					continue;
				}

				EigenType currDiffVec = relevantPoints[i]-relevantPoints[i-1];
				EigenType a,e,d;

				switch (calcCase)
				{
				case 1:
				{
					a = relevantPoints[i-1];
					e = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radInner) * (currDiffVec);
					break;
				}
				case 2:
				{
					a = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radOuter) * (currDiffVec);
					e = relevantPoints[i];
					break;

				}
				case 3:
				{
					a = relevantPoints[i-1];
					e = relevantPoints[i];
					break;
				}
				case 4:
				{
					std::cout << "WARNING: You should not have come so far with calcCase 4\n";
					exit(1);
					break;
				}
				case 5:
				{
					a = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radOuter) * (currDiffVec);
					e = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radInner) * (currDiffVec);
					break;
				}
				default:
				{
					std::cout << "WARNING: No proper calculation case could be determined\n";
					exit(1);
					break;
				}


				} // End of switch

				d = e-a;

				//
				#ifdef DPC_CONT_USE_LOOKUP
					frozenTemp += ((a-this->position).norm()-(e-this->position).norm())* diffRadInv *  look_up_intergral(a,this->position,d);
				#else
					frozenTemp += diffRadInv *  calculate_integral(a,this->position,d);
				#endif

				#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
				if (std::isnan(frozenTemp[0]) || std::isnan(frozenTemp[1]) )
				{
					std::cout << "ERROR: NaN value in full calculation of frozen value, case: " << calcCase << std::endl;
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "outerPoint: " << a.transpose() << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "innerPoint: " << e.transpose() << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "outerDist: " << (a-this->position).norm() << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "innerDist: " << (e-this->position).norm() << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "position: " << this->position.transpose() << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "radInner: " << this->radInner << "\n";
					std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
						<< "radOuter: " << this->radOuter << "\n";
					exit(1);
				}
				#endif


			}// end for loop



			frozenPart += frozenTemp * this->density;



			return;
		}

	// calculate magnetization contributed by given part of the path
	void calculate_partial_value(const EigenType & outerPoint, const EigenType& innerPoint,
			const double & distanceOuter, const double & distanceInner,
			const EigenType & direction,
			EigenType& frozenPart, EigenType& unfrozenPart) const
	{
		if (abs(distanceOuter-distanceInner) < DISTANCE_EPS)
		{
			unfrozenPart += (distanceInner - radInner) * diffRadInv * direction * this->density ;
			return;
		}
		if (distanceOuter <= radInner)
		{
			return;
		}

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (distanceOuter < this->radInner)
		{
			std::cout << "WARNING: Partial value requested for points within inner sphere\n";
			errorCode = 1;
			return;
		}
		if (distanceOuter<distanceInner)
		{
			std::cout << "WARNING: Wrong order of distances for calculation of partial value\n";
			errorCode = 2;
			return;
		}
		#endif


		if (std::abs(distanceOuter-this->radInner) < EQUALITY_EPS)
		{
			return;
		}

		int calcCase = 4;

		if ((distanceOuter > radOuter) && (distanceInner < radInner)) // evolution passing surface in one step
		{
			calcCase = 5;

		}
		else if ((distanceOuter > radOuter) && (distanceInner >= radInner)) // partly unfrozen
		{
			calcCase = 2;
			unfrozenPart += (distanceInner - radInner) * diffRadInv * direction * this->density ;
		}
		else if ((distanceOuter <= radOuter) && (distanceInner < radInner)) // evolution started in ring
		{
			calcCase = 1;
		}
		else //if ((distanceOuter <= radOuter) && (distanceInner >= radInner))
		{
			calcCase = 3;
			unfrozenPart += (distanceInner - radInner) * diffRadInv * direction *this->density;
		}


		EigenType currDiffVec = innerPoint-outerPoint;
		EigenType a,e,d;
		switch (calcCase)
		{
		case 1:
		{
			a = outerPoint;
			e = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radInner) * (currDiffVec);
			break;
		}
		case 2:
		{
			a = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radOuter) * (currDiffVec);
			e = innerPoint;
			break;
		}
		case 3:
		{
			a = outerPoint;
			e = innerPoint;
			break;
		}
		case 4:
		{
			std::cout << "WARNING: You should not have come so far with calcCase 4\n";
			break;
		}
		case 5:
		{
			a = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radOuter) * (currDiffVec);
			e = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radInner) * (currDiffVec);
			break;
		}
		default:
		{
			std::cout << "WARNING: No proper calculation case could be determined\n";
			break;
		}


		} // End of switch

		d = e-a;

		#ifdef DPC_CONT_USE_LOOKUP
			frozenPart += ((a-this->position).norm()-(e-this->position).norm())* diffRadInv * look_up_intergral(a,this->position,d) *this->density;
		#else
			frozenPart += diffRadInv *  calculate_integral(a,this->position,d) *this->density;
		#endif

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (std::isnan(frozenPart[0]) || std::isnan(frozenPart[1]) )
		{
			std::cout << "ERROR: NaN value in partial calculation of frozen value, case: " << calcCase << std::endl;
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "outerPoint: " << outerPoint.transpose() << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "innerPoint: " << innerPoint.transpose() << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "outerDist: " << distanceOuter << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "innerDist: " << distanceInner << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "position: " << this->position.transpose() << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "radInner: " << this->radInner << "\n";
			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< "radOuter: " << this->radOuter << "\n";
			errorCode = 3;
			return;
		}

		errorCode = 0;
		#endif

	}

	// value if path consists of one point
	void calculate_partial_value_single_point(const EigenType& direction,
	const double& distance, EigenType& frozenPart, EigenType& unfrozenPart) const
		{
			if ((distance >radInner))
			{
				// EigenType midVec = -(radOuter-distances[indexOuter]) * diffRadInv * this->position.normalized();
				unfrozenPart += (distance - radInner) * diffRadInv * this->density * direction;
				frozenPart += calculate_zero_state(distance) * this->density;
			}
			else
			{
				frozenPart += calculate_zero_state(radInner)* this->density;
			}

			#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
			if (std::isnan(frozenPart[0]) || std::isnan(frozenPart[1]) )
			{
				std::cout << "ERROR: NaN value in partial calculation of frozen value, initialization case, single point" << std::endl;
				exit(1);
			}
			#endif
			return;
		}



	// calculate initialization value

#ifndef CS_SPHERE_CONT_USE_IDEALIZED_INITIALIZATION
// OLD INITIALIZATION USING x-axis
	EigenType calculate_zero_state(const double& outerDist) const
	{
		EigenType testDir;
		bool zeroCase = false;

		// determine case: y neg, y pos, y zero or x and y zero
		if (this->position[0] < -INIT_DISTANCE_EPS)
		{
			testDir = {-1,0};
		}
		else if (this->position[0] > INIT_DISTANCE_EPS)
		{
			testDir = {1,0};
		}
		else
		{
			if (std::abs(this->position[1]) < INIT_DISTANCE_EPS)
			{
				return EigenType::Zero();
			}
			else
			{
				testDir = {1,0};
				zeroCase = true;
			}
		}

		// construct relevant path
		EigenType startPoint = {this->position[0],0};

		EigenType a = startPoint + testDir *calculate_lambda(startPoint,this->position,testDir, this->radOuter);
		EigenType e = startPoint + testDir *calculate_lambda(startPoint,this->position,testDir, outerDist);

		EigenType d = e-a;

		//
		#ifdef DPC_CONT_USE_LOOKUP
		EigenType res = ((a-this->position).norm()-(e-this->position).norm())* diffRadInv * look_up_intergral(a,this->position,d);
		#else
		EigenType res =  diffRadInv *  calculate_integral(a,this->position,d);
		#endif

		// hysteron origin is on the y-axis => x-component of initial state vector is zero
		if (zeroCase)
		{
			res[0] = 0;
		}

		// Debugging checks
		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (std::isnan(res[0]) || std::isnan(res[1]) )
		{
			std::cout << "ERROR: NaN value in partial calculation of frozen value, initialization case" << std::endl;
			exit(1);
		}
		#endif

		// retrun result
		//std::cout << "norm result:" << res.norm() << std::endl;
		return res;

	}
#else
	EigenType calculate_zero_state(const double& outerDist) const
	{
		if ((std::abs(this->position[0]) < EQUALITY_EPS) && (std::abs(this->position[1]) < EQUALITY_EPS))
		{
			return EigenType::Zero();
		}

		EigenType testDir = EigenType::Zero()-this->position;
		testDir = testDir.normalized();
		EigenType a = testDir * this->radOuter;
		EigenType e = testDir * std::max(this->radInner, outerDist);
		EigenType d = e-a;

		EigenType res =  diffRadInv *  calculate_integral(a,this->position,d);
		return res;


	}
#endif
	// functions to evaluate the magnetization integral (look up works but has no real advantage)
 	EigenType calculate_integral(const EigenType& a, const EigenType& m, const EigenType& d) const
 	{
 		EigenType am = a-m;
 		double amx = am[0];
 		double amy = am[1];
 		double dx = d[0];
 		double dy = d[1];

 		EigenType em = am +d;

 		double dsqInv = d[0]*d[0] + d[1]*d[1];

 		double eNorm = em.norm();
 		double aNorm = am.norm();

 		if (dsqInv < dpcCont::EQUALITY_EPS)
 		{
 			return EigenType::Zero();
 		}

 		dsqInv = 1/dsqInv;

 		double p = 2 *(amx*dx + amy*dy) * dsqInv;
 		double q = (amx*amx + amy*amy) * dsqInv;
 		double root = 4*q-p*p;


 		double Ax = amx*dy*dy - dx*dy*amy;
 		double Ay = amy*dx*dx - dx*dy*amx;

 		double logVal =log(eNorm/aNorm);//;  //eNorm/aNorm -1;//log_taylor_series<1>(eNorm/aNorm);log(eNorm/aNorm);//

 		if (root < EQUALITY_EPS)
 		{
 			EigenType result = {-dx - dsqInv * Ax *  logVal,
 								-dy - dsqInv * Ay *  logVal};
 			return result;
 		}

 		double facx = -dx*std::abs(amx*dy-amy*dx);
 		double facy = -dy*std::abs(amx*dy-amy*dx);

 		double atanVal = hystLib::functions::acos_x( am.dot(em) / (aNorm*eNorm));

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
 		if (std::isnan(atanVal)|| std::isnan(logVal) || std::isinf(logVal))
 		{
 			std::cout << "ERROR: NaN value in intergral calculation, logVal: " << logVal << ", acosVal: " << atanVal << "\n";
 			exit(1);
 		}
		#endif

 		EigenType result ={	-dx - (dsqInv)*(Ax * logVal + (facx) * atanVal),
 							-dy - (dsqInv)*(Ay * logVal + (facy) * atanVal)};

 		return result;
 	}

 	EigenType look_up_intergral (const EigenType& a, const EigenType& m, const EigenType& d) const
 	{
 		if (d.norm() < dpcCont::EQUALITY_EPS)
 		{
 			return EigenType::Zero();
 		}

 		EigenType result;
 		EigenType am = a-m;

 		double refVal = 1/am.norm();
 		am = am * refVal;
 		EigenJacob rotMat;
 		rotMat << am[0], am[1], -am[1], am[0];
 		EigenType dm = rotMat *d * refVal;

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
 		if (dm.norm() >1)
 		{
 			std::cout << "ERROR: Look-up in continuous sphere won't word \n";
 			exit(1);
 		}
		#endif

 		if(dm[1]<0)
 		{
 			// int index0 = int(-dm[0] * dpcCont::conversionMultiplier);
 			// int index1 = int(-dm[1]*dpcCont::conversionMultiplier);
 			// index =  index0 * indexMult + index1;

 			int index = int(-dm[0] * dpcCont::conversionMultiplier) * dpcCont::indexMultiplier
 							+ int(-dm[1]*dpcCont::conversionMultiplier);

 			/*
 			double coeff[8];
 			double sum = 0;
 			for (int i = 0; i < 8; i++)
 			{
 				coeff[i] = dpcCont::contIntTable[index][i];
 				sum += coeff[i];
 			}

			if (std::abs(sum) < 1e-12)
			{
				std::cout << "Fishy coefficients" << index << "\n";
			}
			*/
 			result[0] = dpcCont::contIntTable[index][0] * dm[0] - dm[1] *dpcCont::contIntTable[index][1] -
 					dm[0]*dm[1]*dpcCont::contIntTable[index][2] + dpcCont::contIntTable[index][3];

 			result[1] = dpcCont::contIntTable[index][4] * dm[0] + dm[1] * dpcCont::contIntTable[index][5] +
 					dm[0]*dm[1]*dpcCont::contIntTable[index][6] + dpcCont::contIntTable[index][7];

 		}
 		else
 		{

 			//int index0 = int(-dm[0]*100);
 			//int index1 = int( dm[1]*100);

 			int index = int(-dm[0] * dpcCont::conversionMultiplier) * dpcCont::indexMultiplier
						+ int(dm[1]*dpcCont::conversionMultiplier);

 			result[0] = dpcCont::contIntTable[index][0] * dm[0] + dm[1] *dpcCont::contIntTable[index][1]
					  + dm[0]*dm[1]*dpcCont::contIntTable[index][2] + dpcCont::contIntTable[index][3];

 			result[1] = dpcCont::contIntTable[index][4] * dm[0] + dm[1] *dpcCont::contIntTable[index][5]
						  + dm[0]*dm[1]*dpcCont::contIntTable[index][6] + dpcCont::contIntTable[index][7];
 		}

 		rotMat << am[0], -am[1], am[1], am[0];

 		return rotMat * result;

 	}

 	double calculate_lambda(const EigenType& a, const EigenType& m, const EigenType& d, const double& r ) const
 	{
 		double c1 = d[0]*d[0] + d[1]*d[1];
 		double c2 = 2*(d[0]*(a[0]-m[0]) + d[1]*(a[1]-m[1]));
 		double c3 = (a[0]-m[0])*(a[0]-m[0])+ (a[1]-m[1])*(a[1]-m[1]) - r*r;
 		double root =  c2*c2-4*c3*c1;
 		if (root < 0)
 		{
			#ifdef ENABLE_CS_SPHERE_CONT_SQRT_WARNING
 			std::cout << "WARNING: negative root in lambda calculation in CS_sphere_cont, root was: " << root << std::endl;
			#endif
 			return (-c2 / (2*c1));
 		}
 		else
 		{
 			return (-c2 - std::sqrt(root))/(2*c1);
 		}
 	}


	// global Jacobian calculation functions (unfreezing version should be equivalent)
	EigenJacob calculate_freezing_jacobian(const EigenType& outerPoint, const EigenType& innerPoint,
			const double& distanceOuter, const double& distanceInner) const
	{
		EigenJacob Jacobian = EigenJacob::Zero();

		bool shiftingStart = false;
		bool shiftingEnd = false;


		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (distanceOuter < this->radInner)
		{
			std::cout << "WARNING: Freezing Jacobian requested for segment within the sphere\n";
			return EigenJacob::Zero();
		}
		#endif

		if ((abs(distanceOuter-distanceInner) < EQUALITY_EPS) || (distanceOuter<distanceInner) || (distanceOuter< radInner))
		{
			//std::cout << "WARNING: In calculate Jacobian... <^>-----(*_*)-----<^> \n";
			return Jacobian;
		}

		int calcCase = 4;

		if ((distanceOuter > radOuter) && (distanceInner < radInner)) // evolution passing surface in one step
		{
			shiftingEnd = true;
			shiftingStart = true;
			calcCase = 5;
		}
		else if ((distanceOuter > radOuter) && (distanceInner >= radInner)) // partly unfrozen
		{
			shiftingStart = true;
			calcCase = 2;
			EigenType dist = innerPoint - this->position;
			EigenJacob JacPart2;
			double drdx = dist[0]/dist.norm();
			double drdy = dist[1]/dist.norm();
			JacPart2 << dist.normalized()[0]  *drdx, dist.normalized()[0]  *drdy, dist.normalized()[1]  *drdx, dist.normalized()[1]  *drdy;
			Jacobian +=  diffRadInv *JacPart2;
			Jacobian += (distanceInner - radInner) * diffRadInv * this->calculate_unit_vector_jacobian(innerPoint)/this->density;
		}
		else if ((distanceOuter <= radOuter) && (distanceInner < radInner)) // evolution started in ring
		{
			shiftingEnd = true;
			calcCase = 1;
		}
		else // if ((distanceOuter <= radOuter) && (distanceInner >= radInner))
		{
			calcCase = 3;
			EigenType dist = innerPoint - this->position;
			EigenJacob JacPart2;
			double drdx = dist[0]/dist.norm();
			double drdy = dist[1]/dist.norm();
			JacPart2 << dist.normalized()[0]  *drdx, dist.normalized()[0]  *drdy, dist.normalized()[1]  *drdx, dist.normalized()[1]  *drdy;
			Jacobian +=  diffRadInv *JacPart2;
			Jacobian += (distanceInner - radInner) * diffRadInv * this->calculate_unit_vector_jacobian(innerPoint)/this->density;
		}


		EigenType currDiffVec = innerPoint-outerPoint;
		EigenType a,e,d;
		switch (calcCase)
		{
		case 1:
		{
			a = outerPoint;
			e = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radInner) * (currDiffVec);
			break;
		}
		case 2:
		{
			a = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radOuter) * (currDiffVec);
			e = innerPoint;
			break;

		}
		case 3:
		{
			a = outerPoint;
			e = innerPoint;
			break;
		}
		case 4:
		{
			std::cout << "WARNING: You should not have come so far with calcCase 4\n";
			break;
		}
		case 5:
		{
			a = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radOuter) * (currDiffVec);
			e = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radInner) * (currDiffVec);
			break;
		}
		default:
		{
			std::cout << "WARNING: No proper calculation case could be determined\n";
			break;
		}


		} // End of switch

		d = e-a;
		double 	dMxdex, dMxdey,
				dMydex, dMydey,

				dMxdax, dMxday,
				dMydax, dMyday,

				dexdx, dexdy,
				deydx, deydy,

				daxdx, daxdy,
				daydx, daydy;

		if (shiftingStart)
		{
			calculate_x_y_factor(daxdx, daxdy, daydx, daydy, this->radOuter, outerPoint, innerPoint);
		}
		else
		{
			daxdx = 0.0;
			daxdy = 0.0;
			daydx = 0.0;
			daydy = 0.0;
		}

		if (shiftingEnd)
		{
			calculate_x_y_factor(dexdx, dexdy, deydx, deydy, this->radInner, outerPoint, innerPoint);
		}
		else
		{
			dexdx = 1.0;
			dexdy = 0.0;
			deydx = 0.0;
			deydy = 1.0;
		}

		calculate_Jacobian_coefficients(dMxdex, dMxdey, dMydex, dMydey,
										dMxdax, dMxday, dMydax, dMyday,
										a, e, d, shiftingStart);

		EigenJacob partJac;

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (	std::isnan(dMxdex) || std::isnan(dMxdey) || std::isnan(dMydex) || std::isnan(dMydey) ||
				std::isnan(dMxdax) || std::isnan(dMxday) || std::isnan(dMydax) || std::isnan(dMyday) ||
				std::isnan(dexdx) || std::isnan(dexdy) || std::isnan(deydx) || std::isnan(deydy) ||
				std::isnan(daxdx) || std::isnan(daxdy) || std::isnan(daydx) || std::isnan(daydy)
				)
		{
			std::cout << "Jacobian coefficient is nan\n";
			exit(1);

		}
		#endif

		partJac <<	dMxdex * dexdx + dMxdey * deydx +  dMxdax * daxdx +  dMxday * daydx,
					dMxdex * dexdy + dMxdey * deydy +  dMxdax * daxdy +  dMxday * daydy,
					dMydex * dexdx + dMydey * deydx +  dMydax * daxdx +  dMyday * daydx,
					dMydex * dexdy + dMydey * deydy +  dMydax * daxdy +  dMyday * daydy;

		Jacobian += diffRadInv * partJac;
		//Jacobian += look_up_jacobian(a, this->position, d, e);
		return Jacobian*this->density;
	}

	EigenJacob calculate_unfreezing_jacobian_simple(const EigenType& frozenPoint,
			const EigenType& unfrozenPoint) const
	{

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (std::abs(this->distance(frozenPoint)-this->distance(unfrozenPoint)) > 1e-6)
		{
			std::cout << "WARNING: directional derivative for points with unequal distance requested:  "
					<< this->distance(unfrozenPoint) << " " << this->distance(frozenPoint) << "\n";
			errorCode = 6;
		}
		#endif

		EigenType unfrozenDirection = unfrozenPoint - this->position;
		unfrozenDirection = unfrozenDirection.normalized();

		EigenType frozenDirection = frozenPoint - this->position;
		frozenDirection = frozenDirection.normalized();

		EigenType diffVec = unfrozenDirection - frozenDirection;
		if (diffVec.norm() > dpcCont::EQUALITY_EPS)
		{
			EigenJacob JacobianDir;

			JacobianDir << 	unfrozenDirection[0]*diffVec[0], unfrozenDirection[1]*diffVec[0],
							unfrozenDirection[0]*diffVec[1], unfrozenDirection[1]*diffVec[1];

			JacobianDir *= this->density * diffRadInv ;

			return JacobianDir;
		}
		else
		{
			return EigenJacob::Zero();
		}
	}

	EigenJacob calculate_unfreezing_jacobian(const EigenType& unfrozenPoint, const EigenType& innerPoint, const EigenType& outerPoint,
			const double& distanceInner, const double& distanceOuter) const
	{
		EigenJacob Jacobian = EigenJacob::Zero();
		//EigenType frozenPart = EigenType::Zero();

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (distanceOuter < this->radInner)
		{
			std::cout << "WARNING: Freezing Jacobian requested for segment within the sphere\n";
			return EigenJacob::Zero();
		}
		#endif

		if ((abs(distanceOuter-distanceInner) < EQUALITY_EPS) || (distanceOuter<distanceInner) || (distanceOuter< radInner))
		{
			//std::cout << "WARNING: In calculate Jacobian... <^>-----(*_*)-----<^> \n";
			return Jacobian;
		}

		EigenType a,e,d, currDiffVec;


		currDiffVec = innerPoint -outerPoint;

		if ((distanceOuter > radOuter) && (distanceInner > radInner)) // partly unfrozen
		{
			a = outerPoint + calculate_lambda(outerPoint,this->position, currDiffVec, this->radOuter) * (currDiffVec);
		}
		else
		{
			a = outerPoint;
		}

		e = innerPoint;
		a = a- this->position;
		e = e -this->position;
		d = e-a;
		double dexdL, deydL, dLdr, drdx, drdy, r;
		dexdL= d[0];
		deydL = d[1];

		r = e.norm();


		double 	dMxdex, dMxdey,	dMydex, dMydey,	dMxdax, dMxday,	dMydax, dMyday;
		dLdr = - r / std::sqrt(std::pow(a.dot(d),2)-a.squaredNorm()*d.squaredNorm() + r*r * d.squaredNorm());

		if (std::isnan(dLdr))
		{
			std::cout << "ERROR: Derivative of r with respect to lambda is nan  sqrt term:" << std::pow(a.dot(d),2)-a.squaredNorm()*d.squaredNorm() + r*r * d.squaredNorm() <<"\n";
			return calculate_directional_derivative(innerPoint,
					unfrozenPoint, EigenType::Zero());
		}

		dMxdax=0.0;
		dMxday=0.0;
		dMydax=0.0;
		dMyday=0.0;
		calculate_Jacobian_coefficients(dMxdex, dMxdey, dMydex, dMydey,
										dMxdax, dMxday, dMydax, dMyday,
										a+this->position, e+this->position, d, false);

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (	std::isnan(dMxdex) || std::isnan(dMxdey) || std::isnan(dMydex) || std::isnan(dMydey) ||
				std::isnan(dMxdax) || std::isnan(dMxday) || std::isnan(dMydax) || std::isnan(dMyday)
				)
		{
			std::cout << "ERROR: Unfreezing Jacobian coefficient is nan\n";
			return calculate_directional_derivative(innerPoint,
					unfrozenPoint, EigenType::Zero(), errorCode);

		}
		#endif

		EigenType dist = unfrozenPoint - this->position;
		drdx = dist[0]/dist.norm();
		drdy = dist[1]/dist.norm();

		EigenJacob couplpingJac ;
		couplpingJac << dexdL*dLdr*drdx,
						dexdL*dLdr*drdy,
						deydL*dLdr*drdx,
						deydL*dLdr*drdy;
		Jacobian << dMxdex, dMxdey, dMydex, dMydey;

		Jacobian = Jacobian*couplpingJac;

		EigenJacob JacPart2;
		JacPart2 << dist.normalized()[0]  *drdx, dist.normalized()[0]  *drdy, dist.normalized()[1]  *drdx, dist.normalized()[1]  *drdy;
		Jacobian +=  JacPart2;

		Jacobian *= this->density*diffRadInv;

		return Jacobian;

	}

	#ifndef CS_SPHERE_CONT_USE_IDEALIZED_INITIALIZATION
 //OLD INITIALIZATION USING x-axis
	EigenJacob calculate_unfreezing_jacobian_initial(const EigenType& unfrozenPoint) const
	{
		// reconstruct frozen point
		int calcCase = 0;

		EigenType frozenPoint;
		double dist  = this->distance(unfrozenPoint);

		EigenType testDir = EigenType::Zero();

		if ((this->position[0] < -INIT_DISTANCE_EPS))
		{
			testDir = {-1,0};
		}
		else if ((this->position[0] > INIT_DISTANCE_EPS) )
		{
			testDir = {1,0};
		}
		else
		{
				if ((std::abs(this->position[1]) < INIT_DISTANCE_EPS) )
				{
					calcCase = 2;
				}
				else
				{
					testDir = {1,0};
					calcCase = 1;
				}
		}

		EigenType startPoint = {this->position[0],0};

		frozenPoint = startPoint + testDir *calculate_lambda(startPoint,this->position,testDir, dist);

		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (std::abs(this->distance(frozenPoint)-this->distance(unfrozenPoint)) > 1e-6)
		{
			std::cout << "WARNING: directional derivative for points with unequal distance requested:  "
					<< this->distance(unfrozenPoint) << " " << this->distance(frozenPoint) << "\n";
			errorCode = 7;
		}
		#endif

		// calculate corresponding directions
		EigenType unfrozenDirection = unfrozenPoint - this->position;
		unfrozenDirection = unfrozenDirection.normalized();

		EigenType frozenDirection = frozenPoint - this->position;
		frozenDirection = frozenDirection.normalized();


		if(calcCase ==1)
		{
			frozenDirection[0] = 0;
		}
		else if (calcCase == 2)
		{
			frozenDirection = EigenType::Zero();
		}

		// calculate Jacobian
		EigenType diffVec = unfrozenDirection - frozenDirection;
		EigenJacob JacobianDir;

		JacobianDir << 	unfrozenDirection[0]*diffVec[0], unfrozenDirection[1]*diffVec[0],
						unfrozenDirection[0]*diffVec[1], unfrozenDirection[1]*diffVec[1];

		JacobianDir *= this->density * diffRadInv;

		return JacobianDir;
	}
#else
EigenJacob calculate_unfreezing_jacobian_initial(const EigenType& unfrozenPoint, std::size_t& errorCode) const
{
	EigenType frozenDirection;
	if ((std::abs(this->position[0]) < EQUALITY_EPS) && (std::abs(this->position[1]) < EQUALITY_EPS))
	{
		frozenDirection = EigenType::Zero();

	}
	else
	{
			frozenDirection = EigenType::Zero()-this->position;
			frozenDirection = frozenDirection.normalized();
	}




	EigenType unfrozenDirection = unfrozenPoint - this->position;
	unfrozenDirection = unfrozenDirection.normalized();

	EigenType diffVec = unfrozenDirection - frozenDirection;

	EigenJacob JacobianDir;

	JacobianDir << 	unfrozenDirection[0]*diffVec[0], unfrozenDirection[1]*diffVec[0],
					unfrozenDirection[0]*diffVec[1], unfrozenDirection[1]*diffVec[1];

	JacobianDir *= this->density * diffRadInv ;

	return JacobianDir;


}
	#endif

	// calculation of factors for Jacobian entries
	void calculate_x_y_factor(double& dpxdex, double& dpxdey, double& dpydex, double& dpydey, const double& r,
			const EigenType& a, const EigenType& e) const

	{
		  EigenType d = e-a;

		  double d_x = d[0];
		  double d_y = d[1];

		  double dNormSq = d.squaredNorm();

		  EigenType aR = a-this->position;

		  double a_x = aR[0];
		  double a_y = aR[1];


		  double aNormSq = aR.squaredNorm();

		  double aDd = aR.dot(d);
		  double aDdSq = aDd * aDd;



		  double pdNSq1 = 1.0 / dNormSq;
		  double pdNSq2 = pdNSq1 * pdNSq1;
		  double pdNSq3 = pdNSq2 * pdNSq1;
		  double sqrtVal2 = sqrt(pdNSq2*aDdSq - (aNormSq-r*r)*pdNSq1);

			#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
			  if ((pdNSq2*aDdSq - (aNormSq-r*r)*pdNSq1) < 0)
			  {
				  std::cout << "ERROR: Requested derivative of intersection coordinates, "
						  "but there is no intersection sqrtVal: " << pdNSq2*aDdSq - (aNormSq-r*r)*pdNSq1 << std::endl;
				  exit(1);
			  }
			#endif

		  double sqrtVal1 = 1.0 / sqrtVal2;

		  dpxdex = d_x*(-a_x*pdNSq1+(sqrtVal1*(pdNSq3*(2.0*d_x)*aDdSq - a_x*pdNSq2*aDd - pdNSq2*d_x*(aNormSq-r*r)))+pdNSq2*(2.0*d_x)*aDd)- aDd*pdNSq1 - sqrtVal2;


		  dpxdey = d_x*(-a_y*pdNSq1+(sqrtVal1*(pdNSq3*(2.0*d_y)*aDdSq - a_y*pdNSq2*aDd - pdNSq2*d_y*(aNormSq-r*r)))+pdNSq2*(2.0*d_y)*aDd);


		  dpydex = d_y*(-a_x*pdNSq1+(sqrtVal1*(pdNSq3*(2.0*d_x)*aDdSq - a_x*pdNSq2*aDd - pdNSq2*d_x*(aNormSq-r*r)))+pdNSq2*(2.0*d_x)*aDd);


		  dpydey = d_y*(-a_y*pdNSq1+(sqrtVal1*(pdNSq3*(2.0*d_y)*aDdSq - a_y*pdNSq2*aDd - pdNSq2*d_y*(aNormSq-r*r)))+pdNSq2*(2.0*d_y)*aDd)-aDd*pdNSq1 - sqrtVal2;


	}

	void calculate_Jacobian_coefficients(	double& dMxdex, double& dMxdey, double& dMydex, double& dMydey,
											double& dMxdax, double& dMxday, double& dMydax, double& dMyday,
											const EigenType& a, const EigenType& e, const EigenType& d, const bool& shiftingBegin) const
	{

		EigenType aN = a - this->position;
		EigenType eN = e - this->position;

		const double a_x = aN[0];
		const double a_y = aN[1];

		const double e_x = eN[0];
		const double e_y = eN[1];

		const double d_x = d[0];
		const double d_y = d[1];

		const double d_xSq = d_x*d_x;
		const double d_ySq = d_y*d_y;


		const double aNorm =aN.norm();
		const double eNorm =eN.norm();

		const double aNormInv = 1/aNorm;
		const double eNormInv = 1/eNorm;

		const double aNormSq =aN.squaredNorm();
		const double eNormSq =eN.squaredNorm();
		const double dNormSq =d.squaredNorm();
		const double dNormSqInv =1/dNormSq;

		const double eDa = eN.dot(aN);
		double aXd = a_x*d_y - a_y*d_x;

		const double logVal = log(eNormSq/aNormSq);
		const double acosVal = hystLib::functions::acos_x(eDa / (aNorm*eNorm));

		const double c2 = 1.0/pow(d_xSq + d_ySq,2.0);
		const double c3 = (a_x*d_ySq - a_y*d_x*d_y);
		const double c4 = (a_y*d_xSq - a_x*d_x*d_y);

		double aXdSgn = (aXd>0)-(aXd<0);
		aXd = std::abs(aXd);



		dMxdex = dNormSqInv * ( 						a_y*d_y*0.5*logVal	- e_x/eNormSq*c3 + aXd*acosVal 	- a_y*d_x*aXdSgn*acosVal - d_x*(a_x - e_x*eDa*eNormInv*eNormInv)) 	+ c2*d_x*(c3*logVal - 2.0*d_x*aXd*acosVal) - 1.0;
		dMxdey = dNormSqInv * (  			 -(a_x*d_y-a_y*d_x*0.5)*logVal	- e_y/eNormSq*c3 				+ a_x*d_x*aXdSgn*acosVal - d_x*(a_y - e_y*eDa*eNormInv*eNormInv)) 	+ c2*d_y*(c3*logVal - 2.0*d_x*aXd*acosVal);
		dMydex = dNormSqInv * (	 			 -(a_y*d_x-a_x*d_y*0.5)*logVal	- e_x/eNormSq*c4				- a_y*d_y*aXdSgn*acosVal - d_y*(a_x - e_x*eDa*eNormInv*eNormInv)) 	+ c2*d_x*(c4*logVal - 2.0*d_y*aXd*acosVal);
		dMydey = dNormSqInv * (							a_x*d_x*0.5*logVal	- e_y/eNormSq*c4 + aXd*acosVal 	+ a_x*d_y*aXdSgn*acosVal - d_y*(a_y - e_y*eDa*eNormInv*eNormInv)) 	+ c2*d_y*(c4*logVal - 2.0*d_y*aXd*acosVal) - 1.0;

		if (shiftingBegin)
		{
		dMxdax = dNormSqInv * ( 				-0.5*(a_y*d_y+d_ySq)*logVal	+ a_x/aNormSq*c3 - aXd*acosVal 	+ e_y*d_x*aXdSgn*acosVal - d_x*(e_x - a_x*eDa*aNormInv*aNormInv)) 	- c2*d_x*(logVal*c3 - 2.0*acosVal*d_x*aXd) + 1.0;
		dMxday = dNormSqInv * ((d_x*d_y*0.5 - a_y*d_x*0.5 + a_x*d_y)*logVal + a_y/aNormSq*c3 				- e_x*d_x*aXdSgn*acosVal - d_x*(e_y - a_y*eDa*aNormInv*aNormInv))   - c2*d_y*(logVal*c3 - 2.0*acosVal*d_x*aXd);
		dMydax = dNormSqInv * ((d_x*d_y*0.5 - a_x*d_y*0.5 + a_y*d_x)*logVal + a_x/aNormSq*c4 				+ e_y*d_y*aXdSgn*acosVal - d_y*(e_x - a_x*eDa*aNormInv*aNormInv)) 	- c2*d_x*(logVal*c4 - 2.0*acosVal*d_y*aXd);
		dMyday = dNormSqInv * (					-0.5*(a_x*d_x+d_xSq)*logVal + a_y/aNormSq*c4 - aXd*acosVal 	- e_x*d_y*aXdSgn*acosVal - d_y*(e_y - a_y*eDa*aNormInv*aNormInv)) 	- c2*d_y*(logVal*c4 - 2.0*acosVal*d_y*aXd) + 1.0;
		}

	}

 	// never worked properly
 	EigenJacob look_up_jacobian (const EigenType& a, const EigenType& m, const EigenType& d, const EigenType& e) const
 	{
 			EigenJacob jacobian;
 			EigenType result;
 	 		EigenType am = a-m;
 	 		EigenType em = (a+d)-m;
 	 		double refVal = 1/am.norm();
 	 		am = am * refVal;
 	 		em = em *refVal;
 	 		EigenJacob rotMat;
 	 		rotMat << am[0], am[1], -am[1], am[0];
 	 		EigenType dm = rotMat *d * refVal;
 	 		int index;

 			#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
 	 		if (dm.norm() >1)
 	 		{
 	 			std::cout << "ERROR: Look-up in continuous sphere won't work \n";
 	 			exit(1);
 	 		}
 			#endif

 	 		if(dm[1]<0)
 	 		{

 	 			index = int(-dm[0] * dpcCont::conversionMultiplier) * dpcCont::indexMultiplier
 	 							+ int(-dm[1]*dpcCont::conversionMultiplier);


 	 			jacobian << dpcCont::contIntTable[index][0]  - dm[1]*dpcCont::contIntTable[index][2],
						   -dpcCont::contIntTable[index][1]  - dm[0]*dpcCont::contIntTable[index][2],
							dpcCont::contIntTable[index][4]  + dm[1]*dpcCont::contIntTable[index][6],
							dpcCont::contIntTable[index][5]  + dm[0]*dpcCont::contIntTable[index][6];

 	 		}
 	 		else
 	 		{


 	 			index = int(-dm[0] * dpcCont::conversionMultiplier) * dpcCont::indexMultiplier
 							+ int(dm[1]*dpcCont::conversionMultiplier);

 	 			jacobian << dpcCont::contIntTable[index][0]  + dm[1]*dpcCont::contIntTable[index][2],
							dpcCont::contIntTable[index][1]  + dm[0]*dpcCont::contIntTable[index][2],
							dpcCont::contIntTable[index][4]  + dm[1]*dpcCont::contIntTable[index][6],
							dpcCont::contIntTable[index][5]  + dm[0]*dpcCont::contIntTable[index][6];
 	 		}

 	 		double edx = -1/em.norm() * em[0];
 	 		double edy = -1/em.norm() * em[1];

 	 		EigenType M =  diffRadInv * look_up_intergral(a,this->position,d);

 	 		EigenJacob jacPart;
 	 		jacPart << edx*M[0], edy*M[0], edx*M[1], edy*M[0];

 	 		return ((a-this->position).norm()-(e-this->position).norm())* diffRadInv * rotMat.transpose() * jacobian * rotMat
 	 				+ rotMat.transpose() * jacPart * rotMat;

 	}

 	// legacy implementation
 	void calculate_frozen_part_leg(const std::vector<EigenType> & relevantPoints, const std::vector<EigenType>& directions,
			const std::vector<double>& distances, const std::vector<double>& angles, const EigenType& Hnew, EigenType& result) const
	{


 		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (distances[distances.size()-1] > radOuter)
		{
			std::cout << "WARNING: Surface called as frozen, but it is not, radOuter is "
					<< radOuter << " , queried distance is " << distances.back()
					<< " path point: " << relevantPoints.back().transpose()
					<< " position: " << this->position.transpose()<< "\n";
			result += directions.back() * this->density;
			return;
		}
		if ( (relevantPoints.size() != directions.size()) || (relevantPoints.size() != distances.size()) )
		{
			std::cout << "WARNING: Different size of given distances, points and directions\n";
		}
		#endif

		// determine innermost and outermost relevant point
		const int sizeOfVectors = distances.size()-1;
		int indexOuter = sizeOfVectors;
		int indexInner = sizeOfVectors;
		for (int i = sizeOfVectors; i >= 0; --i)
		{
			if (distances[i] < radInner)
			{
				indexInner = i;
			}

			if (distances[i] > radOuter)
			{
				indexOuter = i;
				break;
			}
			else
			{
				indexOuter = i;
			}
		}


		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (indexOuter > indexInner)
		{
			std::cout << "WARNING: Are you sure about the indices?\n";
		}
		#endif

		double unfrozenPart = 0.0;

		// if the hysteron has never been fully unfrozen, assume demagnetized state (sum over all hysterons)
		if (indexOuter == indexInner)
		{
			if ((distances[indexInner] >radInner))
			{
				unfrozenPart = (distances[indexInner] - radInner) * diffRadInv ;
				result += (directions.back() * unfrozenPart + calculate_zero_state(distances[indexInner]))* this->density;
			}
			else
			{
				result += calculate_zero_state(radInner)* this->density;
			}
			return;
		}

		EigenType addToResult = EigenType::Zero();
		int calcCase = 4;

		for (int i = indexInner; i > indexOuter; i--)
		{
			calcCase = 0;
			if ( (i == indexInner) && ((i-1) != indexOuter))
			{
				if (distances[i] < radInner) // evolution ends outside of the ring
				{
					calcCase = 1;
				}
				else	// segment within ring
				{
					calcCase = 3;
					unfrozenPart = (distances[indexInner] - radInner) * diffRadInv;
					addToResult += unfrozenPart * directions.back();
				}
			}
			else if((i == indexInner) && ((i-1) == indexOuter))
			{
				if ((distances[indexOuter] >= radOuter) && (distances[indexInner] < radInner)) // evolution passing surface in one step
				{
					calcCase = 5;

				}
				else if ((distances[indexOuter] >= radOuter) && (distances[indexInner] > radInner)) // partly unfrozen
				{
					calcCase = 2;
					unfrozenPart = (distances[indexInner] - radInner) * diffRadInv ;
					addToResult += unfrozenPart * directions.back();
				}
				else if ((distances[indexOuter] < radOuter) && (distances[indexInner] < radInner)) // evolution started in ring
				{
					addToResult += calculate_zero_state(distances[indexOuter]);
					//addToResult +=  -(radOuter-distances[indexOuter]) * diffRadInv * this->position.normalized();
					calcCase = 1;
				}
				else
				{
					calcCase = 3;
					unfrozenPart = (distances[indexInner] - radInner) * diffRadInv;
					addToResult += unfrozenPart * directions.back() + calculate_zero_state(distances[indexOuter]);
					//addToResult += unfrozenPart * directions.back() -(radOuter-distances[indexOuter]) * diffRadInv * this->position.normalized();
				}
			}
			else if ((i != indexInner) && ((i-1) == indexOuter))
			{
				if (distances[indexOuter] > radOuter) // evolution passing outer surface
				{
					calcCase = 2;
				}
				else //partly unfrozen
				{

					addToResult += calculate_zero_state(distances[indexOuter]);
					//addToResult +=  -(radOuter-distances[indexOuter]) * diffRadInv * this->position.normalized();
					calcCase = 3;
				}
			}
			else // segment is within the ring
			{
				calcCase = 3;
			}

			if (abs(distances[i-1] -distances[i]) < EQUALITY_EPS)
			{
				continue;
			}

			EigenType currDiffVec = relevantPoints[i]-relevantPoints[i-1];
			EigenType a,e,d;

			switch (calcCase)
			{
			case 1:
			{
				a = relevantPoints[i-1];
				e = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radInner) * (currDiffVec);
				break;
			}
			case 2:
			{
				a = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radOuter) * (currDiffVec);
				e = relevantPoints[i];
				break;

			}
			case 3:
			{
				a = relevantPoints[i-1];
				e = relevantPoints[i];
				break;
			}
			case 4:
			{
				std::cout << "WARNING: You should not have come so far with calcCase 4\n";
				break;
			}
			case 5:
			{
				a = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radOuter) * (currDiffVec);
				e = relevantPoints[i-1] + calculate_lambda(relevantPoints[i-1],this->position, currDiffVec, this->radInner) * (currDiffVec);
				break;
			}
			default:
			{
				std::cout << "WARNING: No proper calculation case could be determined\n";
				break;
			}


			} // End of switch

			d = e-a;

			addToResult += diffRadInv *  calculate_integral(a,this->position,d);

			#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
			if (std::isnan(addToResult[0]) ||std::isnan(addToResult[1]))
			{
				std::cout << "NaN in end result of CS_sphere_cont, last case: " << calcCase << " \n";
				exit(1);
			}
			#endif


		}// end for loop



		#ifdef CS_SPERE_CONT_ENABLE_DEBUGGING_CHECKS
		if (addToResult.norm() > (1+EQUALITY_EPS))
		{
			//std::cout << "WARNING: The added direction is > 0, norm is " << addToResult.norm() << " case was " << calcCase << " term " << term << std::endl;
			std::cout << "WARNING: The added direction is > 1, norm is " << addToResult.norm() << " case was " << calcCase <<  std::endl;

			exit(1);
		}
		#endif

		result += addToResult * this->density;

		return;

	}
};




#endif /* HYSTERESIS_OPERATOR_CRITICAL_SURFACE_SHAPES_CS_SPHERE_CONT_H_ */
