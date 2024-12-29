/*
 * dpc_grid_contructor.h
 *
 *  Created on: 03.07.2017
 *      Author: Stephan Willerich
 *      The class defines the basic interface for all grid constructors.
 *
 *
 */

#ifndef HEADERS_HYSTERESIS_OPERATOR_DPC_GRID_CONTRUCTOR_H_
#define HEADERS_HYSTERESIS_OPERATOR_DPC_GRID_CONTRUCTOR_H_

#include "Hysteresis_Operator/critical_surface/generic_critical_surface.h"
#include <memory>
#include <Eigen/StdVector>



template <int dimension, typename CS =  generic_critical_surface<dimension>>
class dpc_grid_constructor
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

public:
	using EigenType =  typename CS::EigenType;
	using AllocEigenType = Eigen::aligned_allocator<EigenType>;
	struct point_information {
			const EigenType coord;
			const unsigned int surfIndBegin, surfIndEnd;
			const std::vector<double> partialDensSums;
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

			point_information(EigenType coordIn, int surfBegIn, int surfEndIn, std::vector<double> partialDenSumsIn)
			: coord(coordIn), surfIndBegin(surfBegIn), surfIndEnd(surfEndIn), partialDensSums(partialDenSumsIn)
			{}
		};




	std::shared_ptr<std::vector<std::shared_ptr<const CS>>> surfaces =
				std::make_shared<std::vector<std::shared_ptr<const CS>>>();

	std::vector<point_information> memoryAllocationGrid;
	std::vector<EigenType, AllocEigenType> gridExtension;

protected:
	double gridExt;
	double noInsert;
	double linPart;
	std::string extensionMode;

public:
	inline virtual ~dpc_grid_constructor(){}

	virtual double get_grid_extend() const = 0;

	virtual double get_no_insert() const = 0;

	virtual double get_lin_part() const = 0;

	virtual std::string get_extension_mode() const = 0;

};



#endif /* HEADERS_HYSTERESIS_OPERATOR_DPC_GRID_CONTRUCTOR_H_ */
