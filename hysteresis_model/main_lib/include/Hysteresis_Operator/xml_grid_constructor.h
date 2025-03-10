// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#ifndef HEADERS_HYSTERESIS_OPERATOR_XML_GRID_CONSTRUCTOR_H_
#define HEADERS_HYSTERESIS_OPERATOR_XML_GRID_CONSTRUCTOR_H_

#include "Hysteresis_Operator/dpc_grid_contructor.h"
#include "third_party/tinyxml2.h"
#include <fstream>

template <int dimension, typename CS = generic_critical_surface<dimension>>
class xml_grid_constructor: public dpc_grid_constructor<dimension, CS>
{
public:
	using EigenType = typename generic_critical_surface<dimension>::EigenType ;
	typedef typename dpc_grid_constructor<dimension, CS>::point_information point_information;

protected:
	std::stringstream reader;
	unsigned int surfaceCounter = 0;



public:
	xml_grid_constructor(std::string filename);

	double get_grid_extend() const;

	double get_no_insert() const;

	double get_lin_part() const;

	std::string get_extension_mode() const;


protected:
	void read_extension_data (tinyxml2::XMLHandle& intPoint);

	void read_point_data(tinyxml2::XMLHandle& memPoint);

	void clear_reader();

};




#endif /* HEADERS_HYSTERESIS_OPERATOR_XML_GRID_CONSTRUCTOR_H_ */
