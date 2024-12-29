/*
 * continuous_dpc_model.h
 *
 *  Created on: 28.11.2018
 *      Author: Stephan Willerich
 */

#ifndef CONTINUOUS_DPC_MODEL_H_
#define CONTINUOUS_DPC_MODEL_H_

#define ENABLE_DPC_VISUALIZATION_FUNCTIONS
//#define ENABLE_DPC_SQRT_WARNING
//#define ENALBE_DPC_CONT_TIMING
#define ENABLE_DPC_IDENTIFCATION_FUNCTIONS
//#define ENABLE_DPC_DEBUGGING_FUNCTIONS
#define ENABLE_DPC_STATISTICS_FUNCTIONS

#include "Hysteresis_Operator/critical_surface/shapes/CS_sphere_cont.h"
#include "Hysteresis_Operator/malloc_grid.h"
#include "Hysteresis_Operator/Misc_Tools/hyst_config.h"
#include <Eigen/LU>

#ifdef ENABLE_DPC_VISUALIZATION_FUNCTIONS
#include <fstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include "third_party/tinyxml2.h"
#include "hystLib_constants.h"
#endif

#ifdef ENALBE_DPC_CONT_TIMING
#include "cxx_misc_tools/Timer.h"
extern Timer dpcContTimer;
#endif

#include <memory>
#include <Eigen/StdVector>

#include <fstream>



template <int dimension>
class continuous_dpc_model
{
public:
	using sphereType = CS_sphere_cont<dimension>;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)
	using EigenType =  Eigen::Matrix<double, dimension, 1>;
	using EigenJacob = Eigen::Matrix<double, dimension, dimension>;
	using EigenStdVector = std::vector<EigenType,Eigen::aligned_allocator<EigenType>>;


protected:
	// pointer to hysterons
	std::vector<std::shared_ptr<const CS_sphere_cont<dimension>>> contSpheres;

	// pointer to vector with surfaces to calculate values in outside the hysteron area
	const std::shared_ptr< const std::vector<unfrozen_surface<dimension>>> combinedSurfaces;

	// Variable storing the history of H
	EigenStdVector Hevol;

public:


	// malloc grid defining the model
	const std::shared_ptr<const malloc_grid<dimension, CS_sphere_cont<dimension>>> mallocGrid;

	// actual magnetic values
	EigenType Hnew;

	EigenType B;

	// Jacobian of B with respect to H
	EigenJacob JacBH;


	int numberOfHysterons = 0;



private:
	// small structure to store information about inserted points
	struct insertPoint
	{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)
		std::size_t start;
		std::size_t end;
		EigenType point;
		insertPoint (std::size_t startIn, std::size_t endIn, EigenType pointIn)
			: start(startIn), end(endIn), point(pointIn)
		{

		}
		insertPoint (const insertPoint& insertPointIn)
		: start(insertPointIn.start), end(insertPointIn.end), point(insertPointIn.point)
		{

		}
		insertPoint& operator = (const insertPoint& insertPointIn)
		{

			start = insertPointIn.start;
			end = insertPointIn.end;
			point = insertPointIn.point;
			return *this;
		}
	};

	// Variables representing the memory of the hysterons
	EigenStdVector curVals;
	EigenStdVector oldVals;

	std::vector<bool> oldFrozen;
	std::vector<bool> frozen;

	std::size_t timeSteps;

	std::vector<std::size_t> globalExtIndices;
	std::vector<std::size_t> oldGlobalExtIndices;

	std::vector<std::vector<insertPoint>> insertPoints;
	std::vector<std::vector<insertPoint>> oldInsertPoints;





	#ifdef ENALBE_DPC_CONT_TIMING
		double TimeForPath = 0.0;
		double TimeForCalc = 0.0;
		double TimeForPathVec = 0.0;
	#endif

	#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
		int inversionCounter = 0;
		int lineSearchCounter = 0;
		int lineSearchIterationCounter = 0;
		int failedInversions = 0;
		int failedLineSearch = 0;
		double maxDeviation = 0.0;
	#endif
public:

	// constructors
	continuous_dpc_model(const continuous_dpc_model<dimension>& modelIn)
	: contSpheres(modelIn.contSpheres),
	  combinedSurfaces(modelIn.combinedSurfaces),
	  Hevol(modelIn.Hevol),
	  mallocGrid(modelIn.mallocGrid),
	  Hnew(modelIn.Hnew),
	  B(modelIn.B),
	  JacBH(modelIn.JacBH),

	  numberOfHysterons(modelIn.numberOfHysterons),

	  curVals(modelIn.curVals),
	  oldVals(modelIn.oldVals),

	  oldFrozen(modelIn.oldFrozen),
	  frozen(modelIn.frozen),

	  timeSteps(modelIn.timeSteps),

	  globalExtIndices(modelIn.globalExtIndices),
	  oldGlobalExtIndices(modelIn.oldGlobalExtIndices),

	  insertPoints(modelIn.insertPoints),
	  oldInsertPoints(modelIn.oldInsertPoints)
	{

	}

	continuous_dpc_model(const std::shared_ptr<const malloc_grid<dimension, CS_sphere_cont<dimension>>>& malloc_gridIN)
	: combinedSurfaces(malloc_gridIN->get_pointers_to_saturation_surfaces()), mallocGrid(malloc_gridIN), Hnew(EigenType::Zero()),  B(EigenType::Zero())
	{
		contSpheres = initialize_cont_spheres();
		initialize_model();
	}

	// getter functions
	EigenType get_H() const
	{
		return Hnew;
	}

	EigenType get_M() const
	{
		return B*hystLib::constants::mu0inv - Hnew;
	}

	EigenType get_B() const
	{
		return B;
	}

	EigenType get_J() const
	{
		return B-hystLib::constants::mu0* Hnew;
	}

	EigenJacob approximate_jacobian_H_B() const
	{
		return JacBH.inverse();
	}

	EigenJacob approximate_jacobian_B_H() const
	{
		return JacBH;
	}

	int get_number_of_hysterons()const
	{
		return contSpheres.size();
	}

	// reset function
	void reset_all_hysterons()
	{
		Hevol.clear();
		frozen.clear();
		oldFrozen.clear();
		curVals.clear();
		oldVals.clear();
		insertPoints.clear();
		globalExtIndices.clear();
		initialize_model();
	}

	// initialization
	void initialize_model()
	{
		//insertPoints.reserve(combinedSurfaces->size());
		Hevol.push_back(EigenType::Zero());
		Hevol.push_back(EigenType::Zero());
		timeSteps = 1;
		for (std::size_t i = 0; i < combinedSurfaces->size(); i++)
		{
			globalExtIndices.push_back(0);
			insertPoints.push_back(std::vector<insertPoint>());
		}

		oldGlobalExtIndices = globalExtIndices;
		oldInsertPoints = std::vector<std::vector<insertPoint>>(insertPoints);

		for (std::size_t i=0; i<contSpheres.size(); i++)
		{
			frozen.push_back(false);
			oldFrozen.push_back(false);
			curVals.push_back(EigenType::Zero());
		}


		std::vector<std::size_t> freezing;
		std::vector<std::size_t> unfreezing;
		update_insert_points(freezing, unfreezing);

		std::vector<double> distanceEvol;
		distanceEvol.push_back(0.0);

		for (std::size_t i = 0; i < mallocGrid->pointCollection.size(); i++)
		{
			double distance;
			EigenType direction;
			(*combinedSurfaces)[i].distance_and_direction(Hevol.back(),direction, distance);

			int start = mallocGrid->pointCollection[i].surfIndBegin;
			int end = mallocGrid->pointCollection[i].surfIndEnd;

			distanceEvol[0] = distance;



			for (int j = end; j >= start; j--)
			{
				if (distance > contSpheres[j]->radOuter)
				{
					for (int k = j; k>= start; k--) // reset all affected frozen states
					{
						curVals[k] = EigenType::Zero();
						frozen[k] = false;
					}

					break;

				}
				frozen[j] = ((distance) < contSpheres[j]->radInner);

				EigenType tempVal;
				contSpheres[j]->calculate_partial_value_single_point(direction,
						distance, curVals[j], tempVal);
			}
		}

		EigenType test =  EigenType::Zero();
		for(std::size_t i=0; i<contSpheres.size(); i++)
		{
			test = test + curVals[i];
		}

		oldFrozen = frozen;
		oldVals = curVals;


		// JacBH is not initialized yet, there fore an (arbitrary) inverse calculation is performed
		// TODO: Less dirty initialization

		EigenType  Bold = {0.001,0};
		calculate_evolution_B(Bold);
		//approximate_jacobian_B_H_explicitly();


	}

private:

	// history management
	void update_H_evolution(const EigenType& HnewIn)
	{
		this->Hnew = HnewIn;
		Hevol.back() = Hnew;
	}

public:
	void accept_evolution()
	{
		oldGlobalExtIndices = globalExtIndices;
		oldInsertPoints = std::vector<std::vector<insertPoint>>(insertPoints);
		Hevol.push_back(Hnew);
		oldFrozen = frozen;
		oldVals = curVals;
		timeSteps += 1;
	}

	// calculation functions
	void calculate_evolution_H (const EigenType& HnewIn)
	{
		#ifdef ENALBE_DPC_CONT_TIMING
		double tstart = dpcContTimer.elapsed();
		#endif
		if ((HnewIn-Hevol[Hevol.size()-1]).norm() < 1e-13)
		{
			return;
		}
		update_H_evolution(HnewIn);

		std::vector<std::size_t> freezing;
		std::vector<std::size_t> unfreezing;
		update_insert_points(freezing, unfreezing);

		#ifdef ENALBE_DPC_CONT_TIMING
		TimeForPath += dpcContTimer.elapsed() - tstart;
		#endif


		calculate_evolution(HnewIn, freezing, unfreezing);
	}

	void calculate_evolution_B(const EigenType& Bnew)
	{
		#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
		if (std::isnan(Bnew[0])||std::isnan(Bnew[1]))
		{
			std::cout << "ERROR: Requested Bnew is NaN\n";
			exit(1);
		}
		#endif

		EigenType deltaB = Bnew-B;
		double normDeltaB = deltaB.norm();

		if (normDeltaB < dpcCont::TOLERANCE_INVERSION)
		{
			return;
		}

		// push in the "right" direction 1e-6*mu0 = 1.2e-12
		EigenType deltaH = 1e-6 * deltaB / normDeltaB;


		Hnew = Hnew + deltaH;

		std::vector<std::size_t> freezing;
		std::vector<std::size_t> unfreezing;


		// perform first step, pathMatrix and freezing then possess all necessery information to proceed

		#ifdef ENALBE_DPC_CONT_TIMING
		double tstart = dpcContTimer.elapsed();
		#endif
		update_H_evolution(Hnew);
		update_insert_points(freezing, unfreezing);
		#ifdef ENALBE_DPC_CONT_TIMING
		TimeForPath += dpcContTimer.elapsed() - tstart;
		#endif

		calculate_evolution(Hnew, freezing, unfreezing);

		/*
		outStream.open ("results/cont/Data/iterB.txt", std::fstream::app);
		outStream  << std::setprecision(16) << B.transpose() << "\n";
		outStream.close();
		 */

		deltaB = Bnew-B;
		normDeltaB = deltaB.norm();

		int iter = 0;


		EigenJacob linearJac;
		//std::size_t errorCode = 0;

		linearJac << mallocGrid->muRp1mu0, 0,0,mallocGrid->muRp1mu0;

		while ((iter < 100))
		{
			iter +=1;
			// calculate Jacobian


			#ifdef ENALBE_DPC_CONT_TIMING
			tstart = dpcContTimer.elapsed();
			#endif
			JacBH = linearJac;



			EigenType direction;
			double distance;


			for (std::size_t i = 0; i < unfreezing.size(); i++)
			{

				(*combinedSurfaces)[unfreezing[i]].distance_and_direction(Hevol.back(),direction, distance);

				EigenJacob UVJ =
						(*combinedSurfaces)[unfreezing[i]].calculate_unit_vector_jacobian(Hnew) /(*combinedSurfaces)[unfreezing[i]].density;

				int start = mallocGrid->pointCollection[unfreezing[i]].surfIndBegin;
				int end = mallocGrid->pointCollection[unfreezing[i]].surfIndEnd;


				for (int j = end; j >= start; j--)
				{

					// unfrozen => unit vector jacobian
					if (distance > contSpheres[j]->radOuter)
					{
						JacBH += ( ((*mallocGrid).pointCollection[unfreezing[i]]).partialDensSums[j - start]) *mallocGrid->Bsat * UVJ;
						break;
					}
					else
					{
						// completely frozen => no contribution to Jacobian
						frozen[j] = ((distance) < contSpheres[j]->radInner);
						if (oldFrozen[j] && frozen[j])
						{
							 JacBH += EigenJacob::Zero(); // for the sake of completeness
						}
						else
						{
							if (!frozen[j])
							{
								if (insertPoints[unfreezing[i]].size() > 0)
								{
									JacBH += contSpheres[j]->calculate_unfreezing_jacobian_simple(
											insertPoints[unfreezing[i]].back().point, Hevol.back()) * this->mallocGrid->Bsat;
									/*
									double distanceOuter =contSpheres[j]->distance(Hevol[insertPoints[unfreezing[i]].back().start]);
									JacBH += contSpheres[j]->calculate_unfreezing_jacobian(Hevol.back(), insertPoints[unfreezing[i]].back().point,
										Hevol[insertPoints[unfreezing[i]].back().start],
										distance, distanceOuter, errorCode)* this->mallocGrid->Bsat;
									*/
									#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
									if (errorCode!=0)
									{
										std::cout << "ERROR: Code " << errorCode <<" at pos. " << 7
												<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
										write_path_visualization_file("bug.xml","debug");
										exit(1);
									}
									#endif
								}
								else
								{
									JacBH += contSpheres[j]->calculate_unfreezing_jacobian_initial(Hevol.back())
											* this->mallocGrid->Bsat;
									#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
									if (errorCode!=0)
									{
										std::cout << "ERROR: Code " << errorCode <<" at pos. " << 8
												<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
										write_path_visualization_file("bug.xml","debug");
										exit(1);
									}
									#endif
								}
								JacBH += contSpheres[j]->density * UVJ
											* (distance-contSpheres[j]->radInner) * contSpheres[j]->diffRadInv *mallocGrid->Bsat;
							}
						}
					}
				}

			} // end of for loop over all unfreezing points;

			for (std::size_t i = 0; i < freezing.size(); i++)
			{

				(*combinedSurfaces)[freezing[i]].distance_and_direction(Hevol.back(),direction, distance);

				int start = mallocGrid->pointCollection[freezing[i]].surfIndBegin;
				int end = mallocGrid->pointCollection[freezing[i]].surfIndEnd;


				for (int j = end; j >= start; j--)
				{

					// unfrozen => unit vector jacobian
					if (distance > contSpheres[j]->radOuter)
					{
						JacBH += ( ((*mallocGrid).pointCollection[freezing[i]]).partialDensSums[j - start]) /(*combinedSurfaces)[freezing[i]].density  *
								(*combinedSurfaces)[freezing[i]].calculate_unit_vector_jacobian(Hnew) *mallocGrid->Bsat;
						break;

					}
					else
					{

						// completely frozen => no contribution to Jacobian
						frozen[j] = ((distance) < contSpheres[j]->radInner);
						if (oldFrozen[j] && frozen[j])
						{
							 JacBH += EigenJacob::Zero(); // for the sake of completeness
						}
						else
						{
							double distanceInner = distance;
							double distanceOuter =(*combinedSurfaces)[freezing[i]].distance(Hevol[Hevol.size()-2]);
							JacBH += mallocGrid->Bsat * contSpheres[j]->calculate_freezing_jacobian(
										Hevol[Hevol.size()-2], Hevol[Hevol.size()-1],
										distanceOuter, distanceInner);
						}

					}
				}

			} // end of for loop over all freezing points;

			// JacBH should now be consistent

			#ifdef ENALBE_DPC_CONT_TIMING
			TimeForCalc += dpcContTimer.elapsed() - tstart;
			#endif

			if (normDeltaB< dpcCont::TOLERANCE_INVERSION)
			{
				break;
			}


			deltaH = JacBH.partialPivLu().solve(deltaB);

			#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
			if (std::isnan(deltaH[0])||std::isnan(deltaH[1]))
			{
				std::cout << "ERROR: Calculated deltaH is NaN\n";
				exit(1);
			}
			#endif

			if (iter >50)
			{
				Hnew = golden_section_search(deltaH, Bnew);
			}
			else
			{
				if ((deltaH.norm()>100))
				{
					deltaH = 100*deltaH.normalized();
				}
				Hnew =  Hnew + deltaH;
			}
/*
			outStream.open ("results/cont/Data/iterH.txt", std::fstream::app);
			outStream  << std::setprecision(16) << Hnew.transpose() << "\n";
			outStream.close();
*/

			#ifdef ENALBE_DPC_CONT_TIMING
			tstart = dpcContTimer.elapsed();
			#endif

			freezing.clear();
			unfreezing.clear();
			update_H_evolution(Hnew);
			update_insert_points(freezing, unfreezing);

			#ifdef ENALBE_DPC_CONT_TIMING
			TimeForPath += dpcContTimer.elapsed() - tstart;
			#endif

			calculate_evolution(Hnew, freezing, unfreezing);

/*
			outStream.open ("results/cont/Data/iterB.txt", std::fstream::app);
			outStream  << std::setprecision(16) << B.transpose() << "\n";
			outStream.close();
*/
			deltaB = Bnew-B;
			normDeltaB = deltaB.norm();
		}

		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
		inversionCounter += iter;
		#endif

		if (normDeltaB> dpcCont::TOLERANCE_INVERSION)
		{
			#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
				failedInversions += 1;
				maxDeviation = std::max(maxDeviation, normDeltaB);
			#endif
			#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
				std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "\n" << Hnew.transpose() << "\n";
				std::cout << "WARNING: Inversion not successful for norm(B): " << Bnew.norm() << ", reached tol.: " << normDeltaB << "\n";
				std::cout << "Breq: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << Bnew.transpose() << std::endl;
				write_inversion_debugging_file("failedInversion.txt", "debug", Bnew);
			#endif
		}
	}

private:
	void calculate_evolution(const EigenType& HnewIn, const std::vector<std::size_t>& freezing, const std::vector<std::size_t>& unfreezing)
	{
		curVals = oldVals;

		#ifdef ENALBE_DPC_CONT_TIMING
		double tstart2 = dpcContTimer.elapsed();
		double tstart = dpcContTimer.elapsed();
		#endif
		#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
		if (freezing.size()+unfreezing.size() != combinedSurfaces->size())
		{
			std::cout << "ERROR: Evolution without state; Entries state vector: " << freezing.size()+unfreezing.size()
					<< " hysteron groups: " << combinedSurfaces->size()  << std::endl;
			exit(1);
		}
		#endif

		EigenType unfrozenPart = EigenType::Zero();
		EigenType frozenPart =  EigenType::Zero();

		EigenStdVector pathVector;
		pathVector.reserve(Hevol.size()+1);
		std::size_t relSize = 0;

		std::vector<double> distanceEvol;
		distanceEvol.reserve(Hevol.size());

		//std::size_t errorCode = 0;

		for (std::size_t i = 0; i < unfreezing.size(); i++)
		{

			EigenType direction;
			double distance;
			(*combinedSurfaces)[unfreezing[i]].distance_and_direction(Hevol.back(),direction, distance);
			// double distanceSq = distance*distance;

			int start = mallocGrid->pointCollection[unfreezing[i]].surfIndBegin;
			int end = mallocGrid->pointCollection[unfreezing[i]].surfIndEnd;

			for (int j = end; j >= start; j--)
			{

				if (distance >= contSpheres[j]->radOuter)
				{
					unfrozenPart += ( ((*mallocGrid).pointCollection[unfreezing[i]]).partialDensSums[j - start]) * (direction);
					for (int k = j; k>= start; k--) // reset all affected frozen states
					{
						curVals[k] = EigenType::Zero();
						frozen[k] = false;
					}
					break;
				}
				else
				{
					frozen[j] = ((distance) < contSpheres[j]->radInner);

					if (oldFrozen[j] && frozen[j])
					{
						frozenPart += curVals[j];
						continue;
					}

					if (frozen[j])
					{
						contSpheres[j]->calculate_partial_value(
								Hevol[insertPoints[unfreezing[i]].back().start],
								insertPoints[unfreezing[i]].back().point,
								contSpheres[j]->distance(Hevol[insertPoints[unfreezing[i]].back().start]), distance, direction, curVals[j], unfrozenPart);
						#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
						if (errorCode!=0)
						{
							std::cout << "ERROR: Code " << errorCode <<" at pos. " << 1
									<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
							write_path_visualization_file("bug.xml","debug");
							exit(1);
						}
						#endif

						frozenPart +=  curVals[j];
					}
					else
					{

						if (insertPoints[unfreezing[i]].size() == 0)
						{
							curVals[j] = EigenType::Zero();
							contSpheres[j]->calculate_partial_value_single_point(direction,
									distance, curVals[j], unfrozenPart);
							frozenPart += curVals[j];
						}
						else
						{


							//from here the current point is assured to be within the ring
							EigenType removeVal = EigenType::Zero();
							EigenType dummy = EigenType::Zero();

							if (oldInsertPoints[unfreezing[i]].size() > 0)
							{
								std::size_t wipeTo = oldInsertPoints[unfreezing[i]].size();
								for (std::size_t k = oldInsertPoints[unfreezing[i]].size()-1; k != std::size_t(-1); k--)
								{
									if (oldInsertPoints[unfreezing[i]][k].end<=insertPoints[unfreezing[i]].back().start)
									{
										break;
									}
									--wipeTo;
									if (oldInsertPoints[unfreezing[i]][k].start == insertPoints[unfreezing[i]].back().start)
									{
										break;
									}
								}


								// if none of the old memory points are wiped out or if all are wiped out, the reference is the new insert point
								std::size_t startIndex = (wipeTo==0||wipeTo == oldInsertPoints[unfreezing[i]].size()) ?
										insertPoints[unfreezing[i]].back().start :
										std::min(oldInsertPoints[unfreezing[i]][wipeTo].start,insertPoints[unfreezing[i]].back().start);

								for (std::size_t k = startIndex; k < Hevol.size()-2; k++)
								{
									double distHSq  = contSpheres[j]->distanceSq(Hevol[k]);
									if (distHSq < contSpheres[j]->radInnerSq)
									{
										break;
									}

									bool skip = false;

									for (std::size_t l = wipeTo; l < oldInsertPoints[unfreezing[i]].size(); l++)
									{

										if(oldInsertPoints[unfreezing[i]][l].start == k)
										{
											if (contSpheres[j]->distanceSq(Hevol[oldInsertPoints[unfreezing[i]][l].start]) >  contSpheres[j]->radInnerSq)
											{
												contSpheres[j]->calculate_partial_value(
																				Hevol[oldInsertPoints[unfreezing[i]][l].start],
																				oldInsertPoints[unfreezing[i]][l].point,
																				contSpheres[j]->distance(Hevol[oldInsertPoints[unfreezing[i]][l].start]),
																				contSpheres[j]->distance(oldInsertPoints[unfreezing[i]][l].point),
																				direction, removeVal, dummy);
												#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
												if (errorCode!=0)
												{
													std::cout << "ERROR: Code " << errorCode <<" at pos. " << 2
															<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
													write_path_visualization_file("bug.xml","debug");
													exit(1);
												}
												#endif

											}
											k = oldInsertPoints[unfreezing[i]][l].end-1;
											wipeTo++;
											skip = true;
										}
									}
									if(skip)
									{
										continue;
									}
									else
									{
										if (distHSq >  contSpheres[j]->radInnerSq)
										{
											contSpheres[j]->calculate_partial_value(Hevol[k],Hevol[k+1],
													std::sqrt(distHSq),contSpheres[j]->distance(Hevol[k+1]),
													direction, removeVal, dummy);
											#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
											if (errorCode!=0)
											{
												std::cout << "ERROR: Code " << errorCode
														<<" at pos. " << 3 << " received from hysteron" << j << " at point " << unfreezing[i] << std::endl;
												write_path_visualization_file("bug.xml","debug");
												exit(1);
											}
											#endif
										}
									}
								}



							}
							else
							{
								for(std::size_t k = insertPoints[unfreezing[i]].back().start; k<Hevol.size()-2; k++)
								{
									double distHSq  = contSpheres[j]->distanceSq(Hevol[k]);
									if (distHSq >  contSpheres[j]->radInnerSq)
									{
										contSpheres[j]->calculate_partial_value(Hevol[k],Hevol[k+1],
												contSpheres[j]->distance(Hevol[k]),contSpheres[j]->distance(Hevol[k+1]),
												direction, removeVal, dummy);
										#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
										if (errorCode!=0)
										{
											std::cout << "ERROR: Code " << errorCode <<" at pos. " << 4
													<< " received from hysteron " << j  << " at point " << unfreezing[i] << std::endl;
											write_path_visualization_file("bug.xml","debug");
											exit(1);
										}
										#endif
									}
								}
							}
							//EigenType temp= EigenType::Zero();

							curVals[j] -= removeVal;
							contSpheres[j]->calculate_partial_value(
									Hevol[insertPoints[unfreezing[i]].back().start],
									insertPoints[unfreezing[i]].back().point,
									contSpheres[j]->distance(Hevol[insertPoints[unfreezing[i]].back().start]),
									distance, direction, curVals[j], unfrozenPart);
							#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
							if (errorCode!=0)
							{
								std::cout << "ERROR: Code " << errorCode <<" at pos. " << 5
										<< " received from hysteron " << j  << " at point " << unfreezing[i] << std::endl;
								write_path_visualization_file("bug.xml","debug");
								exit(1);
							}
							#endif

							frozenPart += curVals[j];
						}
					}
				}
			}
		}



		for (std::size_t i = 0; i < freezing.size(); i++)
		{

			relSize = Hevol.size();

			EigenType direction;
			double distance;
			(*combinedSurfaces)[freezing[i]].distance_and_direction(Hevol.back(),direction, distance);

			int start = mallocGrid->pointCollection[freezing[i]].surfIndBegin;
			int end = mallocGrid->pointCollection[freezing[i]].surfIndEnd;


			for (int j = end; j >= start; j--)
			{

				if (distance > contSpheres[j]->radOuter)
				{
					unfrozenPart += ( ((*mallocGrid).pointCollection[freezing[i]]).partialDensSums[j - start]) * (direction);
					for (int k = j; k>= start; k--) // reset all affected frozen states
					{
						curVals[k] = EigenType::Zero();
						frozen[k] = false;
					}

					break;
				}
				else
				{
					frozen[j] = ((distance) < contSpheres[j]->radInner);
					if ( (oldFrozen[j] && frozen[j]))
					{
						frozenPart += curVals[j];
						continue;
					}
					else
					{

						double dist1 = distance;
						double dist2 = (*combinedSurfaces)[freezing[i]].distance(Hevol[relSize-2]);

						std::size_t ind1, ind2;

						ind1 =relSize-1;
						ind2 =relSize-2;

						contSpheres[j]->calculate_partial_value(Hevol[ind2], Hevol[ind1],
								dist2, dist1, direction, curVals[j], unfrozenPart);
						frozenPart +=  curVals[j];

						#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
						if (errorCode!=0)
						{
							std::cout << "ERROR: Code " << errorCode <<" at pos. " << 6
									<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
							write_path_visualization_file("bug.xml","debug");
							exit(1);
						}
						#endif
					}

				}
			}
			#ifdef ENALBE_DPC_CONT_TIMING
			TimeForCalc += dpcContTimer.elapsed() - tstart;
			#endif

		}


		B = mallocGrid->Bsat * (unfrozenPart + frozenPart) + mallocGrid->muRp1mu0 * HnewIn;
		if (std::isnan(B[0])||std::isnan(B[1]))
		{
			std::cout << "ERROR: Calculated B is NaN\n";
			exit(1);
		}

	}

	// create relevant path
	void update_insert_points(std::vector<bool>& freezing)
	{
		insertPoints = oldInsertPoints;
		//insertPoints = std::vector<std::vector<insertPoint>>(oldInsertPoints);
		globalExtIndices = oldGlobalExtIndices;

		EigenType distVec;
		EigenType diff;
		double lambda;
		double distSq;

		// calculate distance and direction for each point in malloc grid
		for (std::size_t i = 0; i < (*combinedSurfaces).size(); ++i)
		{
			distSq = (Hevol.back()-(*combinedSurfaces)[i].position).squaredNorm();

			// Check if new point is a global maximum with respect to the previous evolution of H
			if( (distSq> (*combinedSurfaces)[i].radMaxSq)
					|| (distSq >= ((*combinedSurfaces)[i].position-Hevol[globalExtIndices[i]]).squaredNorm()) )
			{
				freezing.push_back(false);
				globalExtIndices[i] = timeSteps;
				insertPoints[i].clear();
			}
			else
			{
				std::size_t refIndex = Hevol.size()-2;
				double  refDistSq = (Hevol[refIndex]-(*combinedSurfaces)[i].position).squaredNorm();


				if (distSq-refDistSq > 0) // unfreeze
				{
					freezing.push_back(false);

					// search for last point with greater distance
					// all points after the last inserted points cannot be larger than the current point
					std::size_t startInd;
					if(insertPoints[i].size()<1)
					{
						startInd = Hevol.size()-2;
					}
					else
					{
						if ((insertPoints[i].back().point -(*combinedSurfaces)[i].position).squaredNorm() < distSq)
							startInd = insertPoints[i].back().start;
						else
							startInd = Hevol.size()-2;
					}


					for (std::size_t j = startInd ; j >= globalExtIndices[i] ; j--)
					{
						if (distSq < (Hevol[j]-(*combinedSurfaces)[i].position).squaredNorm())
						{
							refIndex = j;
							break;
						}
					}

					// Calculate inserted point
					diff = Hevol[refIndex+1]- Hevol[refIndex];
					distVec = Hevol[refIndex] - (*combinedSurfaces)[i].position;

					double a = std::pow(diff[0],2) +  std::pow(diff[1],2);
					double b = -2*(diff[0]* distVec[0] + diff[1]* distVec[1]);
					double c = std::pow(distVec[0],2)+ std::pow(distVec[1],2) - distSq;

					lambda = (+b-std::sqrt(b*b - 4*a*c))/(2*a);

					EigenType newPoint = Hevol[refIndex] + lambda *diff;


					// Discard every point that is inserted after the remaining point

					std::size_t j = insertPoints[i].size()-1;
					while ((j!=std::size_t(-1)) && insertPoints[i][j].start >= refIndex)
					{
						j--;
						insertPoints[i].pop_back();
					}

					#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
					if (std::isnan(newPoint[0]) || std::isnan(newPoint[1]))
					{
						std::cout << "ERROR: Point with NaN coordinates insert normal\n";
						exit(1);
					}
					#endif

					insertPoints[i].emplace_back(refIndex, timeSteps, newPoint);

				}
				else if (distSq-refDistSq < 0) // freeze
				{
					diff = Hevol[refIndex+1]- Hevol[refIndex];
					distVec = Hevol[refIndex] - (*combinedSurfaces)[i].position;
					lambda=(std::fmax(-1.0, std::fmin(-distVec.dot(diff)/diff.squaredNorm(), 1.0)));

					// lambda will always be in the interval [-1.0, 1.0], normally lambda should be lambda>=0,
					// however for lambda == 0.0 and lambda == 1.0 no insertion is needed
					if ( lambda > 0.0  && lambda < 1.0)
					{
						lambda += lambda-1.0;
						EigenType newPoint = Hevol[refIndex]  + lambda * diff;
						#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
						if (std::isnan(newPoint[0]) || std::isnan(newPoint[1]))
						{
							std::cout << "ERROR: Point with NaN coordinates in insert extra\n";
							exit(1);
						}
						#endif
						insertPoints[i].emplace_back(timeSteps-1, timeSteps, newPoint);
						freezing.push_back(false);
					}
					else
					{
						freezing.push_back(true);
					}

				}
				else // distance remains (exactly) the same
				{
					freezing.push_back(false);
					// No change...
				}
			}
		}
	} // END of update_insert_points(...)

	void update_insert_points(std::vector<std::size_t>& freezing, std::vector<std::size_t>& unfreezing)
		{
			insertPoints = oldInsertPoints;
			//insertPoints = std::vector<std::vector<insertPoint>>(oldInsertPoints);
			globalExtIndices = oldGlobalExtIndices;

			EigenType distVec;
			EigenType diff;
			double lambda;
			double distSq;

			// calculate distance and direction for each point in malloc grid
			for (std::size_t i = 0; i < (*combinedSurfaces).size(); ++i)
			{
				distSq = (Hevol.back()-(*combinedSurfaces)[i].position).squaredNorm();

				// Check if new point is a global maximum with respect to the previous evolution of H
				if( (distSq> (*combinedSurfaces)[i].radMaxSq)
						|| (distSq >= ((*combinedSurfaces)[i].position-Hevol[globalExtIndices[i]]).squaredNorm()) )
				{
					unfreezing.push_back(i);
					globalExtIndices[i] = timeSteps;
					insertPoints[i].clear();
				}
				else
				{
					std::size_t refIndex = Hevol.size()-2;
					double  refDistSq = (Hevol[refIndex]-(*combinedSurfaces)[i].position).squaredNorm();


					if (distSq-refDistSq >= 0) // unfreeze
					{


						// search for last point with greater distance
						// all points after the last inserted points cannot be larger than the current point
						std::size_t startInd;
						if(insertPoints[i].size()<1)
						{
							startInd = Hevol.size()-2;
						}
						else
						{
							if ((insertPoints[i].back().point -(*combinedSurfaces)[i].position).squaredNorm() < distSq)
								startInd = insertPoints[i].back().start;
							else
								startInd = Hevol.size()-2;
						}

						for (std::size_t j = startInd ; j >= globalExtIndices[i] ; j--)
						{
							if (distSq < (Hevol[j]-(*combinedSurfaces)[i].position).squaredNorm())
							{
								refIndex = j;
								break;
							}
						}

						// Calculate inserted point
						diff = Hevol[refIndex+1]- Hevol[refIndex];
						distVec = Hevol[refIndex] - (*combinedSurfaces)[i].position;

						double a = std::pow(diff[0],2) +  std::pow(diff[1],2);
						double b = -2*(diff[0]* distVec[0] + diff[1]* distVec[1]);
						double c = std::pow(distVec[0],2)+ std::pow(distVec[1],2) - distSq;

						lambda = (+b-std::sqrt(b*b - 4*a*c))/(2*a);

						EigenType newPoint = Hevol[refIndex] + lambda *diff;


						// Discard every point that is inserted after the remaining point

						std::size_t j = insertPoints[i].size()-1;
						while ((j!=std::size_t(-1)) && insertPoints[i][j].start >= refIndex)
						{
							j--;
							insertPoints[i].pop_back();
						}

						#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
						if (std::isnan(newPoint[0]) || std::isnan(newPoint[1]))
						{
							std::cout << "ERROR: Point with NaN coordinates insert normal\n";
							exit(1);
						}
						#endif

						insertPoints[i].emplace_back(refIndex, timeSteps, newPoint);

						unfreezing.push_back(i);

					}
					else if (distSq-refDistSq < 0) // freeze
					{
						diff = Hevol[refIndex+1]- Hevol[refIndex];
						distVec = Hevol[refIndex] - (*combinedSurfaces)[i].position;
						lambda=(std::fmax(-1.0, std::fmin(-distVec.dot(diff)/diff.squaredNorm(), 1.0)));

						// lambda will always be in the interval [-1.0, 1.0], normally lambda should be lambda>=0,
						// however for lambda == 0.0 and lambda == 1.0 no insertion is needed
						if ( lambda > 0.0  && lambda < 1.0)
						{
							lambda += lambda-1.0;
							EigenType newPoint = Hevol[refIndex]  + lambda * diff;
							#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
							if (std::isnan(newPoint[0]) || std::isnan(newPoint[1]))
							{
								std::cout << "ERROR: Point with NaN coordinates in insert extra\n";
								exit(1);
							}
							#endif
							insertPoints[i].emplace_back(timeSteps-1, timeSteps, newPoint);
							unfreezing.push_back(i);
						}
						else
						{
							freezing.push_back(i);
						}

					}
					else // distance remains (exactly) the same
					{
						unfreezing.push_back(i);
						// No change...
					}
				}
			}
		} // END of update_insert_points(...)

	void build_all_relevant_pahts(std::vector<EigenStdVector>& pathMatrix) const
	{
		#ifdef ENALBE_DPC_CONT_TIMING
		double tstart2 = dpcContTimer.elapsed();
		#endif
		for (std::size_t j = 0; j <  (*combinedSurfaces).size(); j++)
		{
			pathMatrix[j].clear();
			if (insertPoints[j].size()==0)
			{

				//pathVector = EigenStdVector(Hevol.begin()+globalExtIndices[indPoint], Hevol.end());

				for (std::size_t i =globalExtIndices[j]; i < Hevol.size(); i++)
				{
					pathMatrix[j].push_back(Hevol[i]);
				}

			}
			else
			{
				std::size_t noOfInserts = insertPoints[j].size()-1;
				std::size_t insertNo = 0;
				std::size_t lookForInd = insertPoints[j][insertNo].start;

				for (std::size_t i =globalExtIndices[j]; i < Hevol.size(); i++)
				{
					pathMatrix[j].push_back(Hevol[i]);
					if (i == lookForInd)
					{
						pathMatrix[j].push_back(insertPoints[j][insertNo].point);

						i = insertPoints[j][insertNo].end-1;

						if (insertNo < noOfInserts)
						{
							insertNo += 1;
							lookForInd = insertPoints[j][insertNo].start;
						}
					}
				}
			}

		}
		#ifdef ENALBE_DPC_CONT_TIMING
		TimeForPathVec += dpcContTimer.elapsed() - tstart2;
		#endif

	}

	void build_relevant_path(EigenStdVector& pathVector, const std::size_t& indPoint) const
	{

		#ifdef ENALBE_DPC_CONT_TIMING
		double tstart2 = dpcContTimer.elapsed();
		#endif
		pathVector.clear();
		if (insertPoints[indPoint].size()==0)
		{
			//pathVector = EigenStdVector(Hevol.begin()+globalExtIndices[indPoint], Hevol.end());

			//pathVector.insert(pathVector.begin(), Hevol.begin()+globalExtIndices[indPoint], Hevol.end());

			for (std::size_t i =globalExtIndices[indPoint]; i < Hevol.size(); i++)
			{
				pathVector.push_back(Hevol[i]);
			}


		}
		else
		{

			std::size_t noOfInserts = insertPoints[indPoint].size()-1;
			std::size_t insertNo = 0;
			std::size_t lookForInd = insertPoints[indPoint][insertNo].start;

			for (std::size_t i =globalExtIndices[indPoint]; i < Hevol.size(); i++)
			{
				pathVector.push_back(Hevol[i]);
				if (i == lookForInd)
				{
					pathVector.push_back(insertPoints[indPoint][insertNo].point);

					i = insertPoints[indPoint][insertNo].end-1;

					if (insertNo < noOfInserts)
					{
						insertNo += 1;
						lookForInd = insertPoints[indPoint][insertNo].start;
					}
				}
			}

		}
		#ifdef ENALBE_DPC_CONT_TIMING
		TimeForPathVec += dpcContTimer.elapsed() - tstart2;
		#endif

	}

	// functions for step-width calculation
	EigenType golden_section_search(const EigenType& deltaH, const EigenType& Bnew)
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
		lineSearchCounter += 1;
		#endif

		double endFac = 1.5;
		double bracketInc = endFac-endFac/ hystLib::constants::goldenRatio;
		double resVec[4];
		double alphaVec[4] = {0,bracketInc, endFac-bracketInc, endFac};

		int iter = 0;

		EigenType Htest;
		const EigenType Horig = Hnew;


		// find bracket
		while ( ((resVec[1]>resVec[0]) || (resVec[1]>resVec[3])) || iter<1 )
		{
			if (iter > 0)
			{
				if (resVec[1]<resVec[0])
				{
					alphaVec[0]=alphaVec[0]+bracketInc;
					alphaVec[1]=alphaVec[1]+bracketInc;
					alphaVec[2]=alphaVec[2]+bracketInc;
					alphaVec[3]=alphaVec[3]+bracketInc;
				}
				else
				{
					endFac = alphaVec[1];
					bracketInc = endFac-endFac / hystLib::constants::goldenRatio;
					alphaVec[0]=0;
					alphaVec[1]=bracketInc;
					alphaVec[2]=endFac-bracketInc;
					alphaVec[3]=endFac;
				}
			}

			iter++;


			Htest = Horig + alphaVec[0] * deltaH;
			calculate_evolution_H(Htest);
			resVec[0] = (B-Bnew).norm();

			Htest = Horig + alphaVec[1] * deltaH;
			calculate_evolution_H(Htest);
			resVec[1] = (B-Bnew).norm();

			Htest = Horig + alphaVec[3] * deltaH;
			calculate_evolution_H(Htest);
			resVec[3] = (B-Bnew).norm();
			if (iter > 10)
			{
				#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
					std::cout << "WARNING: No proper bracket found during line search\n";
				#endif
				#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
					failedLineSearch += 1;
				#endif

				return Horig + deltaH;
			}

		}

		Htest = Horig + alphaVec[2] * deltaH;
		calculate_evolution_H(Htest);
		resVec[2] = (B-Bnew).norm();

		// line search
		iter = 0;
		double tolFac = 1e-1;
		while ((iter < 100) && ((alphaVec[3]-alphaVec[0])>tolFac*(std::abs(alphaVec[1])+ std::abs(alphaVec[2]) )))
		{
			iter ++;
			if (resVec[2] > resVec[1])
			{
				// close in from the right
				alphaVec[3]=alphaVec[2];
				resVec[3]=resVec[2];
				alphaVec[2]=alphaVec[1];
				resVec[2]=resVec[1];

				alphaVec[1]=alphaVec[0]+alphaVec[3]-alphaVec[2];

				Htest = Horig + alphaVec[1] * deltaH;
				calculate_evolution_H(Htest);
				resVec[1] = (B-Bnew).norm();

			}
			else if(resVec[1] > resVec[2])
			{
				// close in from the left
				alphaVec[0]=alphaVec[1];
				resVec[0]=resVec[1];
				alphaVec[1]=alphaVec[2];
				resVec[1]=resVec[2];

				alphaVec[2]=alphaVec[3]-(alphaVec[1]-alphaVec[0]);

				Htest = Horig + alphaVec[2] * deltaH;
				calculate_evolution_H(Htest);
				resVec[2] = (B-Bnew).norm();
			}
		}

		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			lineSearchIterationCounter += iter;
		#endif

		if (resVec[2]>resVec[1])
		{
			return alphaVec[1]*deltaH + Horig;
		}
		else
		{
			return alphaVec[2]*deltaH + Horig;
		}
	}


public:
	double get_max_deviation()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return maxDeviation;
		#else
			return -1;
		#endif
	}

	int get_iteration_number()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return inversionCounter;
		#else
			return -1;
		#endif
	}

	int get_line_search_number()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return lineSearchCounter;
		#else
			return -1;
		#endif
	}

	int get_line_search_iteration_number()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return lineSearchIterationCounter;
		#else
			return -1;
		#endif
	}

	int get_failed_inversions_number()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return failedInversions;
		#else
			return -1;
		#endif
	}

	int get_failed_line_search_number()
	{
		#ifdef ENABLE_DPC_STATISTICS_FUNCTIONS
			return failedLineSearch;
		#else
			return -1;
		#endif
	}
#ifdef ENABLE_DPC_IDENTIFCATION_FUNCTIONS
	void gather_unit_magnetizations(EigenType HnewIn, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result)
	{
		this->Hnew = HnewIn;
		std::vector<bool> freezing;
		update_H_evolution(Hnew);
		update_insert_points(freezing);


		EigenStdVector pathVector;
		for (std::size_t i = 0; i < mallocGrid->pointCollection.size(); i++)
		{
			bool distancesCalculated = false;
			std::vector<double> distanceEvol;
			build_relevant_path(pathVector, i);

			EigenType direction;
			double distance;
			(*combinedSurfaces)[i].distance_and_direction(pathVector.back(),direction, distance);

			/*
			if (distance < dpcCont::EQUALITY_EPS)
			{
				(*combinedSurfaces)[i].distance_and_direction(pathVector[pathVector.size()-2],direction, distance);
			}
*/
			int start = mallocGrid->pointCollection[i].surfIndBegin;
			int end = mallocGrid->pointCollection[i].surfIndEnd;


			for (int j = end; j >= start; j--)
			{
				EigenType unfrozenVal = EigenType::Zero();
				EigenType frozenVal = EigenType::Zero();
				if (distance > contSpheres[j]->radOuter)
				{
					unfrozenVal = (direction);
					frozen[j] = false;
					curVals[j] = EigenType::Zero();
				}
				else
				{
					frozen[j] = ((distance) < contSpheres[j]->radInner);
					if (oldFrozen[j] && frozen[j])
					{
						//std::cout << "Should not be here\n";
						frozenVal += curVals[j];
					}
					else
					{
						if (!freezing[i] || timeSteps == 1)
						{
							if (!distancesCalculated)
							{
								for (std::size_t k = 0; k < pathVector.size(); k++)
								{
									distanceEvol.push_back((*combinedSurfaces)[i].distance(pathVector[k]));
								}
								distancesCalculated = true;
							}
							curVals[j] = EigenType::Zero();
							contSpheres[j]->calculate_frozen_part(pathVector, direction,
									distanceEvol,  Hnew, curVals[j], unfrozenVal);
							frozenVal += curVals[j];
						}
						else
						{


							double dist1 = distance;
							double dist2 = (*combinedSurfaces)[i].distance(pathVector[pathVector.size()-2]);

							std::size_t ind1, ind2;
							if (abs(dist1-dist2) < dpcCont::EQUALITY_EPS)
							{

								ind1 =pathVector.size()-2;
								ind2 =pathVector.size()-3;
								dist1 = dist2;
								dist2 =(*combinedSurfaces)[i].distance(pathVector[ind2]);
							}
							else
							{
								ind1 =pathVector.size()-1;
								ind2 =pathVector.size()-2;
							}
							std::size_t errCode;

							contSpheres[j]->calculate_partial_value(pathVector[ind2], pathVector[ind1],
									dist2, dist1, direction, curVals[j],unfrozenVal);

							frozenVal += curVals[j];
						}

					}

				}
				result.block(0,j,2,1) = unfrozenVal + frozenVal;
			}

		}

		this->accept_evolution();
	}
#endif

	void write_path_visualization_file(const std::string& fileNameRaw, const std::string& folderName) const;

	void write_inversion_debugging_file(const std::string& fileNameRaw, const std::string& folderName, const EigenType& Breq) const
	{
		std::string fileName = (boost::filesystem::path(folderName)/fileNameRaw).string(); 
		//std::string fileName(ewt_tools::concatenate_folders(folderName.c_str(),fileNameRaw.c_str()));
		if (!boost::filesystem::exists(folderName.c_str()))
		{
			boost::filesystem::create_directories(folderName);
		}

		std::ofstream outStream;
		outStream.open (fileName);
		for (std::size_t i = 0; i < Hevol.size()-1; i++)
		{
			outStream  << std::setprecision(16) << Hevol[i].transpose() << "\n";
		}
		outStream << Breq.transpose() << "\n";

		outStream.close();
	}

	EigenJacob approximate_jacobian_B_H_explicitly()
	{
		std::vector<std::size_t> freezing;
		std::vector<std::size_t> unfreezing;
		update_insert_points(freezing, unfreezing);

		EigenJacob linearJac;
		linearJac << mallocGrid->muRp1mu0, 0,0,mallocGrid->muRp1mu0;

		// calculate Jacobian


		#ifdef ENALBE_DPC_CONT_TIMING
		tstart = dpcContTimer.elapsed();
		#endif
		JacBH = linearJac;



		EigenType direction;
		double distance;
		//std::size_t errorCode = 0;


		for (std::size_t i = 0; i < unfreezing.size(); i++)
		{

			(*combinedSurfaces)[unfreezing[i]].distance_and_direction(Hevol.back(),direction, distance);

			EigenJacob UVJ =
					(*combinedSurfaces)[unfreezing[i]].calculate_unit_vector_jacobian(Hnew) /(*combinedSurfaces)[unfreezing[i]].density;

			int start = mallocGrid->pointCollection[unfreezing[i]].surfIndBegin;
			int end = mallocGrid->pointCollection[unfreezing[i]].surfIndEnd;


			for (int j = end; j >= start; j--)
			{

				// unfrozen => unit vector jacobian
				if (distance > contSpheres[j]->radOuter)
				{
					JacBH += ( ((*mallocGrid).pointCollection[unfreezing[i]]).partialDensSums[j - start]) *mallocGrid->Bsat * UVJ;
					break;
				}
				else
				{
					// completely frozen => no contribution to Jacobian
					frozen[j] = ((distance) < contSpheres[j]->radInner);
					if (oldFrozen[j] && frozen[j])
					{
						 JacBH += EigenJacob::Zero(); // for the sake of completeness
					}
					else
					{
						if (!frozen[j])
						{
							if (insertPoints[unfreezing[i]].size() > 0)
							{
								JacBH += contSpheres[j]->calculate_unfreezing_jacobian_simple(
										insertPoints[unfreezing[i]].back().point, Hevol.back()) * this->mallocGrid->Bsat;
								/*
								double distanceOuter =contSpheres[j]->distance(Hevol[insertPoints[unfreezing[i]].back().start]);
								JacBH += contSpheres[j]->calculate_unfreezing_jacobian(Hevol.back(), insertPoints[unfreezing[i]].back().point,
									Hevol[insertPoints[unfreezing[i]].back().start],
									distance, distanceOuter, errorCode)* this->mallocGrid->Bsat;
								*/
								#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
								if (errorCode!=0)
								{
									std::cout << "ERROR: Code " << errorCode <<" at pos. " << 7
											<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
									write_path_visualization_file("bug.xml","debug");
									exit(1);
								}
								#endif
							}
							else
							{
								JacBH += contSpheres[j]->calculate_unfreezing_jacobian_initial(Hevol.back())
										* this->mallocGrid->Bsat;
								#ifdef ENABLE_DPC_DEBUGGING_FUNCTIONS
								if (errorCode!=0)
								{
									std::cout << "ERROR: Code " << errorCode <<" at pos. " << 8
											<< " received from hysteron " << j << " at point " << unfreezing[i] << std::endl;
									write_path_visualization_file("bug.xml","debug");
									exit(1);
								}
								#endif
							}
							JacBH += contSpheres[j]->density * UVJ
										* (distance-contSpheres[j]->radInner) * contSpheres[j]->diffRadInv *mallocGrid->Bsat;
						}
					}
				}
			}

		} // end of for loop over all unfreezing points;

		for (std::size_t i = 0; i < freezing.size(); i++)
		{

			(*combinedSurfaces)[freezing[i]].distance_and_direction(Hevol.back(),direction, distance);

			int start = mallocGrid->pointCollection[freezing[i]].surfIndBegin;
			int end = mallocGrid->pointCollection[freezing[i]].surfIndEnd;


			for (int j = end; j >= start; j--)
			{

				// unfrozen => unit vector jacobian
				if (distance > contSpheres[j]->radOuter)
				{
					JacBH += ( ((*mallocGrid).pointCollection[freezing[i]]).partialDensSums[j - start]) /(*combinedSurfaces)[freezing[i]].density  *
							(*combinedSurfaces)[freezing[i]].calculate_unit_vector_jacobian(Hnew) *mallocGrid->Bsat;
					break;

				}
				else
				{

					// completely frozen => no contribution to Jacobian
					frozen[j] = ((distance) < contSpheres[j]->radInner);
					if (oldFrozen[j] && frozen[j])
					{
						 JacBH += EigenJacob::Zero(); // for the sake of completeness
					}
					else
					{
						double distanceInner = distance;
						double distanceOuter =(*combinedSurfaces)[freezing[i]].distance(Hevol[Hevol.size()-2]);
						JacBH += mallocGrid->Bsat * contSpheres[j]->calculate_freezing_jacobian(
									Hevol[Hevol.size()-2], Hevol[Hevol.size()-1],
									distanceOuter, distanceInner);
					}

				}
			}

		} // end of for loop over all freezing points;

		// JacBH should now be consistent


	return JacBH;
	}
private:
	std::vector<std::shared_ptr<const CS_sphere_cont<dimension>>> initialize_cont_spheres();
};

#endif /* CONTINUOUS_DPC_MODEL_H_ */
