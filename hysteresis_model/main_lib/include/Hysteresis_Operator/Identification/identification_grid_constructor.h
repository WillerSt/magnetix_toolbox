// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HYSTERESIS_OPERATOR_IDENTIFICATION_IDENTIFICATION_GRID_CONSTRUCTOR_H_
#define HYSTERESIS_OPERATOR_IDENTIFICATION_IDENTIFICATION_GRID_CONSTRUCTOR_H_
#include "Hysteresis_Operator/dpc_grid_contructor.h"
#include "Hysteresis_Operator/critical_surface/generic_critical_surface.h"
#include "Hysteresis_Operator/critical_surface/shapes/CS_sphere_gen.h"
template <int dimension>
class identification_grid_constructor: public dpc_grid_constructor<dimension>
{
private:
	typedef typename dpc_grid_constructor<dimension>::EigenType EigenType;
	typedef typename dpc_grid_constructor<dimension>::point_information point_information;
	const double xMax, rMax, rSpacing, xSpacing;

	const int dim = dimension;
public:
	identification_grid_constructor(double xMaxIn, double rMaxIn, int xNumSteps , int rNumSteps)
	: xMax(xMaxIn), rMax(rMaxIn), rSpacing(rMax/rNumSteps), xSpacing(2*xMax/xNumSteps)
	{
		int surfaceCounter = 0;
		switch (dim)
		{
			case 1:
				for (double x = -xMax; x < (xMax + 1e-6); x += xSpacing)
				{
					for (double r = 0; r < (rMax + 1e-6); x += rSpacing)
					{
						// TODO: Implement for one dimension.
					}
				}
				break;

			case 2:
/*
				for (double x = -xMax; x < (xMax + 1e-6); x += xSpacing)
				{
					for  (double y = -xMax; y < (xMax + 1e-6); y += xSpacing)
					{
						point_information memoryPoint;
						EigenType coord;
						coord << x,y;
						std::vector<double> radii;
						std::vector<double> densities;
						memoryPoint.surfIndBegin = surfaceCounter;
						for (double r = 0; r < (rMax + 1e-6); r += rSpacing)
						{
							this->surfaces->push_back(new CS_sphere_gen<2>(coord, r, 1, surfaceCounter));
							radii.push_back(r);
							densities.push_back(1);
							surfaceCounter++;
						}
						memoryPoint.surfIndEnd = surfaceCounter - 1;

						memoryPoint.coord = coord;
						//memoryPoint.radii = radii;
						//memoryPoint.densities = densities;
						this->memoryAllocationGrid.push_back(memoryPoint);
					}
				}
*/

				for (double x = -xMax; x < (xMax + 1e-6); x += xSpacing)
				{
						EigenType coord;
						coord << x,0;
						std::vector<double> radii;
						std::vector<double> densities;
						std::vector<double> partialDensSums;
						int surfIndBegin = surfaceCounter;
						double densSum = 0.0;
						for (double r = 0; r < (rMax + 1e-6); r += rSpacing)
						{
							this->surfaces->push_back(new CS_sphere_gen<2>(coord, r, 1, surfaceCounter));
							radii.push_back(r);
							densities.push_back(1);
							densSum += 1.0;
							partialDensSums.push_back(densSum);
							surfaceCounter++;

						}
						int surfIndEnd = surfaceCounter - 1;

						//memoryPoint.coord = coord;
						//memoryPoint.radii = radii;
						//memoryPoint.densities = densities;
						this->memoryAllocationGrid.push_back(point_information(coord, surfIndBegin, surfIndEnd, partialDensSums));
				}


				/*
				for (double y = -xMax; y < 0; y += xSpacing)
				{
					point_information memoryPoint;
					EigenType coord;
					coord << 0,y;
					std::vector<double> radii;
					std::vector<double> densities;
					memoryPoint.surfIndBegin = surfaceCounter;
					for (double r = 0; r < (rMax + 1e-6); r += rSpacing)
					{
						this->surfaces->push_back(new CS_sphere_gen<2>(coord, r, 1, surfaceCounter));
						radii.push_back(r);
						densities.push_back(1);
						surfaceCounter++;
					}
					memoryPoint.surfIndEnd = surfaceCounter - 1;

					memoryPoint.coord = coord;
					this->memoryAllocationGrid.push_back(memoryPoint);

				}
				for (double y = xSpacing; y < (xMax+1e-6); y += xSpacing)
				{
					point_information memoryPoint;
					EigenType coord;
					coord << 0,y;
					std::vector<double> radii;
					std::vector<double> densities;
					memoryPoint.surfIndBegin = surfaceCounter;
					for (double r = 0; r < (rMax + 1e-6); r += rSpacing)
					{
						this->surfaces->push_back(new CS_sphere_gen<2>(coord, r, 1, surfaceCounter));
						radii.push_back(r);
						densities.push_back(1);
						surfaceCounter++;
					}
					memoryPoint.surfIndEnd = surfaceCounter - 1;

					memoryPoint.coord = coord;
					//memoryPoint.radii = radii;
					//memoryPoint.densities = densities;
					this->memoryAllocationGrid.push_back(memoryPoint);

				}
*/
				break;

			case 3:
				for (double x = -xMax; x < (xMax + 1e-6); x += xSpacing)
				{
					for (double y = -xMax; y < (xMax + 1e-6); y += xSpacing)
					{
						for (double z = -xMax; z < (xMax + 1e-6); z += xSpacing)
						{
							for (double r = 0; r < (rMax + 1e-6); x += rSpacing)
							{
								// TODO: Implement stuff for three dimensions

							}
						}

					}

				}
				break;
		}




	}
	double get_grid_extend() const
	{
		return 0;
	}

	double get_no_insert() const
	{
		return 0;
	}

	std::string get_extension_mode() const
	{
		return "s";
	}

	double get_lin_part() const
	{
		return 0.0;
	}

};


#endif /* HYSTERESIS_OPERATOR_IDENTIFICATION_IDENTIFICATION_GRID_CONSTRUCTOR_H_ */
