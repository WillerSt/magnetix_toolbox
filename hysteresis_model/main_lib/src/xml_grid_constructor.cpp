// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include "Hysteresis_Operator/xml_grid_constructor.h"
#include "boost/filesystem.hpp"
#include "Hysteresis_Operator/critical_surface/known_shapes.h"
#include <iostream>

template <int dimension, typename CS>
xml_grid_constructor<dimension, CS>::xml_grid_constructor(std::string filename)
{
		const char* conversion = filename.c_str();


		int numberOfMemPoints;
		int numberOfIntPoints;
		int dim;

		// Check if file exists return error otherwise
		if (boost::filesystem::exists(filename ))
		{
			tinyxml2::XMLDocument xmlIn;
			xmlIn.LoadFile(conversion);
			tinyxml2::XMLHandle docHandle(xmlIn);
			tinyxml2::XMLHandle currMemPoint(xmlIn);
			tinyxml2::XMLHandle currIntPoint(xmlIn);


			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("NumberOfPoints");
			reader >> numberOfMemPoints;
			clear_reader();


			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("Dim");
			this->reader >> dim;
			clear_reader();

			if (xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("linFac") != nullptr)
			{
				this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("linFac");
				this->reader >> this->linPart;
				clear_reader();
			}
			else
			{
				std::cout << "INFO: No linear part specified, defaulting to muR = 0\n";
				this->linPart = 0;
			}
			// Check if dimension of file and created constructor are matching
			if (dim!=dimension)
			{
				std::cout << "Error: Dimension of DPC-Data (" << dim
						<<") does not match dimension of constructor (" << dimension << ") \n";
				return;
			}
			//std:: cout << "No. of memory points: " << numberOfMemPoints << "\n\n";
			clear_reader();

			// read noInsert value
			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationInfo")->Attribute("noInsert");
			this->reader >> (this->noInsert);
			clear_reader();

			// read memory allocation grid
			currMemPoint = docHandle.FirstChildElement("DPC_Data").FirstChildElement("MemoryAllocationGrid")
					.FirstChildElement("MemoryPoint");
			this->read_point_data(currMemPoint);

			for (int i = 0; i < numberOfMemPoints-1; i++)
			{
				currMemPoint = currMemPoint.NextSibling();
				this->read_point_data(currMemPoint);
			}

			// read interpolation points

			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("NumberOfPoints");
			this->reader >> numberOfIntPoints;
			clear_reader();
			currIntPoint = docHandle.FirstChildElement("DPC_Data").FirstChildElement("InterpolationGrid").
					FirstChildElement("coordinate");
			read_extension_data(currIntPoint);

			for (int i = 0; i < numberOfIntPoints-1; i++)
			{
				currIntPoint = currIntPoint.NextSibling();
				read_extension_data(currIntPoint);
			}

			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("Extend");
			this->reader >> (this->gridExt);
			clear_reader();
			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("Shape");
			this->reader >> (this->extensionMode);
			clear_reader();
		}
		else
		{
			std::cout << "ERROR: Specified XML-File for the Hysteresis Model does not exist or cannot be accessed!\n \n";
			exit(1);
		}

}



template <int dimension, typename CS>
void xml_grid_constructor<dimension,CS>::read_extension_data (tinyxml2::XMLHandle& intPoint)
{
	std::vector<double> coordVec;
	typename dpc_grid_constructor<dimension>::EigenType coordEigen;
	double value;
	tinyxml2::XMLElement* intCoord = intPoint.ToElement();

	reader << intCoord->GetText();
	while (reader >> value)
	{
		coordVec.push_back(value);
	}
	clear_reader();

	for (int i = 0; i<dimension; i++)
	{
		coordEigen[i]=coordVec[i];
	}
	this->gridExtension.push_back(coordEigen);

}



template <int dimension , typename CS>
void xml_grid_constructor<dimension, CS>::read_point_data(tinyxml2::XMLHandle& memPoint)
{
	std::vector<CS*> surfvec;
	double density, value;
	typename dpc_grid_constructor<dimension>::EigenType coordEigen;
	// Get Handles to child nodes containing the point informations

	// read coordinate
	tinyxml2::XMLElement* coordinate = memPoint.FirstChildElement("coordinate").ToElement();
	std::vector<double> coordVec;
	clear_reader();
	reader << coordinate->GetText();
	while (reader >> value)
	{
		coordVec.push_back(value);
	}
	for (int i = 0; i<dimension; i++)
	{
		coordEigen[i]=coordVec[i];
	}
	clear_reader();

	// read critical surface data

	tinyxml2::XMLHandle currSurfaceHandle = memPoint.FirstChildElement("CriticalSurface");
	tinyxml2::XMLElement* critSurf = currSurfaceHandle.ToElement();




	 int beginTemp = surfaceCounter;

	while (critSurf != nullptr)
	{

		reader << critSurf->Attribute("Density");
		reader >> density;
		clear_reader();

		reader << critSurf->Attribute("ShapeIdentifier");
		char Identifier;
		reader >> Identifier;
		clear_reader();

		if(CS::shapeIdentifier != Identifier)
		{
			std::cout << "Error: shape Identifier of template '" << CS::shapeIdentifier << "' does not match identifier in xml \'" << Identifier <<"\'\n";
			exit(1);
		}

		double radius;

		switch (Identifier)
		{

		case 's':
		{
			reader << critSurf->Attribute("Radius");
			reader >> radius;
			clear_reader();
			this->surfaces->push_back(std::make_shared<const CS>(coordEigen, radius, density, surfaceCounter));
			break;
		}
		case 'S':
		{
			reader << critSurf->Attribute("Radius");
			double outerRad;
			double innerRad;
			bool foundRadii= false;
			reader << critSurf->Attribute("Radius");
			reader >> radius;
			clear_reader();

			if (critSurf->Attribute("OuterRadius")!=nullptr)
			{
				reader << critSurf->Attribute("OuterRadius");
				reader >> outerRad;
				clear_reader();
				foundRadii = true;
			}

			if (critSurf->Attribute("InnerRadius") != nullptr)
			{
				reader << critSurf->Attribute("InnerRadius");
				reader >> innerRad;
				clear_reader();
			}
			else
			{
				foundRadii = false;
			}

			if (foundRadii)
			{
				this->surfaces->push_back(std::make_shared<const CS>(coordEigen, radius,  innerRad, outerRad, density, surfaceCounter));
			}
			else
			{
				std::cout << "WARNING: Continuous Sphere specified, but no outer or inner radius are given";
				this->surfaces->push_back(std::make_shared<const CS>(coordEigen, radius, density, surfaceCounter));
			}
			break;
		}
		default:
		{
			std::cout << "Error: Encountered unknown shape identifier " << Identifier << std::endl;
			exit(1);
		}
		}





		unsigned int plannedID;
		this->reader << critSurf->Attribute("ID");
		this->reader >> plannedID;
		this->clear_reader();
		if ((plannedID != (*(this->surfaces)).back()->get_surfaceID()) || (plannedID != this->surfaces->size() - 1))
		{
			std::cout<< "WARNING: Surface ID " << (*(this->surfaces)).back()->get_surfaceID() << " do not match planned IDs " << plannedID << " (only important for identification) \n";
		}

		surfaceCounter ++;
		currSurfaceHandle = currSurfaceHandle.NextSibling();
		critSurf = currSurfaceHandle.ToElement();
		// critSurf->NextSibling();
	}

	int endTemp = surfaceCounter - 1;

	std::vector<double> partialDensSums;
	double sum = 0;
	for (int i = beginTemp; i <= endTemp; ++i)
	{
		sum +=  (*((*this->surfaces)[i])).density;
		partialDensSums.push_back(sum);
	}

	if (int(partialDensSums.size()) != (endTemp-beginTemp+1))
	{
		std::cout << "WARNING: Partial density sums are not consistent\n";
	}
	point_information pointInfo(coordEigen, beginTemp, endTemp, partialDensSums);


	// store memory point
	if (pointInfo.surfIndEnd < pointInfo.surfIndBegin)
	{
		std::cout << "WARNING: Malloc grid point without any associated surfaces encountered.\n";
	}
	else
	{
		this->memoryAllocationGrid.push_back(pointInfo);
	}
}


// Specialization for generic surface type

template <int dimension>
class xml_grid_constructor<dimension, generic_critical_surface<dimension>>:public dpc_grid_constructor<dimension, generic_critical_surface<dimension>>
{

	using EigenType = typename generic_critical_surface<dimension>::EigenType ;
	using point_information = typename dpc_grid_constructor<dimension, generic_critical_surface<dimension>>::point_information;

public:
	xml_grid_constructor(std::string filename)
{
		const char* conversion = filename.c_str();


		int numberOfMemPoints;
		int numberOfIntPoints;
		int dim;

		// Check if file exists return error otherwise
		if (boost::filesystem::exists(filename ))
		{
			tinyxml2::XMLDocument xmlIn;
			xmlIn.LoadFile(conversion);
			tinyxml2::XMLHandle docHandle(xmlIn);
			tinyxml2::XMLHandle currMemPoint(xmlIn);
			tinyxml2::XMLHandle currIntPoint(xmlIn);


			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("NumberOfPoints");
			reader >> numberOfMemPoints;
			clear_reader();

			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("Dim");
			this->reader >> dim;

			if (xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("linFac") != nullptr)
			{
				this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("MemoryAllocationGrid")->Attribute("linFac");
				this->reader >> this->linPart;
				clear_reader();
			}
			else
			{
				std::cout << "INFO: No linear part specified, defaulting to muR = 0\n";
				this->linPart = 0;
			}

			// Check if dimension of file and created constructor are matching
			if (dim!=dimension)
			{
				std::cout << "Error: Dimension of DPC-Data (" << dim
						<<") does not match dimension of constructor (" << dimension << ") \n";
				return;
			}
			//std:: cout << "No. of memory points: " << numberOfMemPoints << "\n\n";
			clear_reader();

			// read noInsert value
			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationInfo")->Attribute("noInsert");
			this->reader >> (this->noInsert);
			clear_reader();

			// read memory allocation grid
			currMemPoint = docHandle.FirstChildElement("DPC_Data").FirstChildElement("MemoryAllocationGrid")
					.FirstChildElement("MemoryPoint");
			this->read_point_data(currMemPoint);

			for (int i = 0; i < numberOfMemPoints-1; i++)
			{
				currMemPoint = currMemPoint.NextSibling();
				this->read_point_data(currMemPoint);
			}

			// read interpolation points

			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("NumberOfPoints");
			this->reader >> numberOfIntPoints;
			clear_reader();
			currIntPoint = docHandle.FirstChildElement("DPC_Data").FirstChildElement("InterpolationGrid").
					FirstChildElement("coordinate");
			read_extension_data(currIntPoint);

			for (int i = 0; i < numberOfIntPoints-1; i++)
			{
				currIntPoint = currIntPoint.NextSibling();
				read_extension_data(currIntPoint);
			}

			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("Extend");
			this->reader >> (this->gridExt);
			clear_reader();
			this->reader << xmlIn.FirstChildElement("DPC_Data")->FirstChildElement("InterpolationGrid")->Attribute("Shape");
			this->reader >> (this->extensionMode);
			clear_reader();
		}
		else
		{
			std::cout << "ERROR: Specified XML-File for the Hysteresis Model does not exist or cannot be accessed!\n \n";
			exit(1);
		}
}

	double get_grid_extend() const {return this->gridExt;}

	double get_no_insert() const {return this->noInsert;}

	std::string get_extension_mode() const {return this->extensionMode;}

	double get_lin_part() const {return this->linPart;}

protected:
	std::stringstream reader;
	unsigned int surfaceCounter = 0;

private:
	void clear_reader()
	{
		this->reader.str(std::string());
		this->reader.clear();
	}

	void read_extension_data (tinyxml2::XMLHandle& intPoint)
	{
		std::vector<double> coordVec;
		typename dpc_grid_constructor<dimension>::EigenType coordEigen;
		double value;
		tinyxml2::XMLElement* intCoord = intPoint.ToElement();

		reader << intCoord->GetText();
		while (reader >> value)
		{
			coordVec.push_back(value);
		}
		clear_reader();

		for (int i = 0; i<dimension; i++)
		{
			coordEigen[i]=coordVec[i];
		}
		this->gridExtension.push_back(coordEigen);

	}

	void read_point_data(tinyxml2::XMLHandle& memPoint)
	{
		std::vector<generic_critical_surface<dimension>*> surfvec;
		double density, value;
		typename dpc_grid_constructor<dimension>::EigenType coordEigen;
		// Get Handles to child nodes containing the point informations

		// read coordinate
		tinyxml2::XMLElement* coordinate = memPoint.FirstChildElement("coordinate").ToElement();
		std::vector<double> coordVec;
		this->clear_reader();
		this->reader << coordinate->GetText();
		while (this->reader >> value)
		{
			coordVec.push_back(value);
		}
		for (int i = 0; i<dimension; i++)
		{
			coordEigen[i]=coordVec[i];
		}
		this->clear_reader();

		// read critical surface data

		tinyxml2::XMLHandle currSurfaceHandle = memPoint.FirstChildElement("CriticalSurface");
		tinyxml2::XMLElement* critSurf = currSurfaceHandle.ToElement();


		 int beginTemp = surfaceCounter;

		while (critSurf != nullptr)
		{

			this->reader << critSurf->Attribute("Density");
			this->reader >> density;
			this->clear_reader();

			this->reader << critSurf->Attribute("ShapeIdentifier");
			char Identifier;
			this->reader >> Identifier;
			this->clear_reader();

			switch (Identifier)
			{

			case 's':
				this->reader << critSurf->Attribute("Radius");
				double radius;
				this->reader >> radius;
				this->clear_reader();
				this->surfaces->push_back(std::make_shared<CS_sphere_gen<dimension>>(coordEigen, radius, density, this->surfaceCounter));
				break;
			default:
				std::cout << "Error: Encountered unknown shape identifier " << Identifier << std::endl;
				exit(1);
			}


			unsigned int plannedID;
			this->reader << critSurf->Attribute("ID");
			this->reader >> plannedID;
			this->clear_reader();
			if ((plannedID != (*(this->surfaces)).back()->get_surfaceID()) || (plannedID != this->surfaces->size() - 1))
			{
				std::cout<< "WARNING: Surface ID " << (*(this->surfaces)).back()->get_surfaceID() << " do not match planned IDs " << plannedID << " (only important for identification) \n";
			}

			critSurf->NextSibling();
			this->surfaceCounter ++;
			currSurfaceHandle = currSurfaceHandle.NextSibling();
			critSurf = currSurfaceHandle.ToElement();
		}

		int endTemp = surfaceCounter - 1;

		std::vector<double> partialDensSums;
		double sum = 0;
		for (int i = beginTemp; i <= endTemp; ++i)
		{
			sum +=  (*((*this->surfaces)[i])).density;
			partialDensSums.push_back(sum);
		}

		if (int(partialDensSums.size())!= (endTemp-beginTemp+1))
		{
			std::cout << "WARNING: Partial density sums are not consistent\n";
		}
		point_information pointInfo(coordEigen, beginTemp, endTemp, partialDensSums);



		// store memory point
		if (pointInfo.surfIndEnd < pointInfo.surfIndBegin)
		{
			std::cout << "WARNING: Malloc grid point without any associated surfaces encountered.\n";
		}
		else
		{
			this->memoryAllocationGrid.push_back(pointInfo);
		}
	}
};

template <int dimension, typename CS>
void xml_grid_constructor<dimension, CS>::clear_reader()
{
	this->reader.str(std::string());
	this->reader.clear();
}


template<int dimension, typename CS>
double xml_grid_constructor<dimension, CS>::get_grid_extend() const
{
	return this->gridExt;
}

template<int dimension, typename CS>
double xml_grid_constructor<dimension, CS>::get_lin_part() const
{
	return this->linPart;
}

template<int dimension, typename CS>
double xml_grid_constructor<dimension, CS>::get_no_insert() const
{
	return this->noInsert;
}




template<int dimension, typename CS>
std::string xml_grid_constructor<dimension, CS>::get_extension_mode() const
{
	return this->extensionMode;
}

template class xml_grid_constructor<1, generic_critical_surface<1>>;
template class xml_grid_constructor<2, generic_critical_surface<2>>;
template class xml_grid_constructor<3, generic_critical_surface<3>>;

template class xml_grid_constructor<1, CS_sphere<1>>;
template class xml_grid_constructor<2, CS_sphere<2>>;
template class xml_grid_constructor<3, CS_sphere<3>>;

template class xml_grid_constructor<1, CS_sphere_cont<1>>;
template class xml_grid_constructor<2, CS_sphere_cont<2>>;
template class xml_grid_constructor<3, CS_sphere_cont<3>>;


