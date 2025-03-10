// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HEADERS_HYSTERESIS_OPERATOR_MALLOC_GRID_
#define HEADERS_HYSTERESIS_OPERATOR_MALLOC_GRID_
#include <Hysteresis_Operator/critical_surface/special/unfrozen_surface.h>
#include "Hysteresis_Operator/dpc_grid_contructor.h"


template<int dimension, typename CS =  generic_critical_surface<dimension>>
class malloc_grid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(true)

private:
	using point_information = typename dpc_grid_constructor<dimension, CS>::point_information;
public:
	typedef typename Eigen::Matrix<double, dimension, 1> EigenType;
	const std::vector<point_information> pointCollection;
	const std::shared_ptr<std::vector<std::shared_ptr<const CS>>> surfaces;
	const std::shared_ptr<std::vector<unfrozen_surface<dimension>>> saturationSurfaces;

public:
	const double Bsat;
	const unsigned int numberOfSurfaces;
	const double muR;
	const double muRmu0;
	const double muRp1mu0;
public:


	malloc_grid(dpc_grid_constructor<dimension, CS>& GridConstructor, double BsatIn);

	void display_all_surface_infos() const;
	void display_reduced_surface_infos() const;
	void display_point_info(int i) const;
	point_information access_grid_point(int i) const;

	unsigned int get_number_of_points() const;
	// unsigned int get_number_of_radii(unsigned int i) const;

	std::shared_ptr<const CS> get_pointer_to_surface(unsigned int i) const;
	std::shared_ptr<std::vector<unfrozen_surface<dimension>>> get_pointers_to_saturation_surfaces() const;

	/*
	std::pair<EigenType, double> get_surface_info_from_ID(unsigned int surfaceID) const;
	*/

private:
	std::vector<point_information> initialize_malloc_grid(dpc_grid_constructor<dimension, CS>& GridConstructor);
	std::shared_ptr<std::vector<unfrozen_surface<dimension>>> construct_saturation_surfaces();

};




#endif /* HEADERS_HYSTERESIS_OPERATOR_MALLOC_GRID_ */
