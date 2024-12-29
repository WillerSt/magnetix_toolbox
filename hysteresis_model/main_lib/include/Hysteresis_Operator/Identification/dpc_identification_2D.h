/*
 * dpc_identification_2D.h
 *
 *  Created on: 25.01.2018
 *      Author: Stephan Willerich
 */

#ifndef HYSTERESIS_OPERATOR_IDENTIFICATION_DPC_IDENTIFICATION_2D_H_
#define HYSTERESIS_OPERATOR_IDENTIFICATION_DPC_IDENTIFICATION_2D_H_

#include <Hysteresis_Operator/Identification/identification_grid_constructor.h>

#include "Hysteresis_Operator/malloc_grid.h"
//
#include "Hysteresis_Operator/xml_grid_constructor.h"
#include <boost/filesystem.hpp>

#include "Hysteresis_Operator/interpolation_grid.h"

#include <hdf5/serial/H5Cpp.h>
#include <hdf5/serial/H5File.h>

#include <Eigen/Sparse>

#include <fstream>
#include <iostream>

template<class dpc_model>
class dpc_idenftification_2D
{
private:
	using EigenType = dpc_grid_constructor<2>::EigenType;

	const double xMax;
	const double rMax;

	const int stepsX;
	const int stepsR;

	const std::string filename;

	int numberOfDataPoints;


	std::vector<std::vector<std::pair<EigenType, EigenType>>> MeasurementSequences;

	std::vector<std::size_t> sequenceSelection;

	std::unique_ptr<dpc_grid_constructor<2, typename dpc_model::sphereType>> gridConst;

	const std::shared_ptr<const malloc_grid<2, typename dpc_model::sphereType>> identGrid
			= std::make_shared<const malloc_grid<2, typename dpc_model::sphereType>>(*gridConst, 1);

	dpc_model identModel = dpc_model(identGrid);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b;


public:
	dpc_idenftification_2D(std::string filenameIN, std::vector<std::size_t> sequenceSelectionIN, const std::string& matrixPath)
		:xMax(1000), rMax(300), stepsX(10), stepsR(50), filename(filenameIN),
		 sequenceSelection(sequenceSelectionIN), gridConst(new identification_grid_constructor<2>(xMax, rMax, stepsX, stepsR))
	{
		load_identification_data();
		construct_matrices();
		write_matrices_to_hdf5(matrixPath);
	};
	dpc_idenftification_2D(std::string filenameIN, double xMaxIn, double rMaxIn, int stepsXIn, int stepsRIn,
			std::vector<std::size_t> sequenceSelectionIN, const std::string& matrixPath)
		:xMax(xMaxIn), rMax(rMaxIn), stepsX(stepsXIn), stepsR(stepsRIn), filename(filenameIN),
		 sequenceSelection(sequenceSelectionIN), gridConst(new identification_grid_constructor<2>(xMax, rMax, stepsX, stepsR))
	{
		load_identification_data();
		construct_matrices();
		write_matrices_to_hdf5(matrixPath);

	};
	dpc_idenftification_2D(std::string filenameIN, std::string gridXMLIn,
			std::vector<std::size_t> sequenceSelectionIN, const std::string& matrixPath)
		:xMax(0), rMax(0), stepsX(0), stepsR(0), filename(filenameIN), sequenceSelection(sequenceSelectionIN),
		 gridConst(new xml_grid_constructor<2, typename dpc_model::sphereType>(gridXMLIn))
	{
		load_identification_data();
		construct_matrices();
		write_matrices_to_hdf5(matrixPath);
	};

private:
	void load_identification_data()
	{

		if (!boost::filesystem::exists(filename))
		{
			std::cout << "ERROR (dpc_identification_2D.h): Specified input file " << filename << " has not been found.\n";
			exit(1);
		}

		int dimIn;

		std::stringstream reader;
		tinyxml2::XMLDocument xmlIn;
		xmlIn.LoadFile(filename.c_str());
		tinyxml2::XMLHandle docHandle(xmlIn);

		tinyxml2::XMLHandle sequence(xmlIn);

		reader << xmlIn.FirstChildElement("HBMeasurement")->Attribute("dim");

		reader >> dimIn;

		if (dimIn != 2)
		{
			std::cout << "WARNING: The data you provided are " << dimIn << "-dimensional, 2-dimensional data are required.\n";
		}

		reader.str(std::string());
		reader.clear();

		reader << xmlIn.FirstChildElement("HBMeasurement")->Attribute("dataPoints");

		reader >> numberOfDataPoints;

		sequence = xmlIn.FirstChildElement("HBMeasurement")->FirstChildElement("Sequence");


		int sequenceNo = 0;

		// TODO Is seuence != nullptr correct?
		while (sequence.ToElement() != nullptr)
		{
			std::vector<std::pair<EigenType, EigenType>> sequenceVals;;

			// std::cout << "Sequence No. " << sequenceNo << std::endl;
			std::stringstream BxText;
			BxText << sequence.FirstChildElement("Bx").ToElement()->GetText();
			std::stringstream ByText;
			ByText << sequence.FirstChildElement("By").ToElement()->GetText();
			std::stringstream HxText;
			HxText << sequence.FirstChildElement("Hx").ToElement()->GetText();
			std::stringstream HyText;
			HyText << sequence.FirstChildElement("Hy").ToElement()->GetText();

			double BxVal, ByVal, HxVal, HyVal;

			while (BxText >>BxVal)
			{
				ByText >> ByVal;
				HxText >> HxVal;
				HyText >> HyVal;

				EigenType B = {BxVal, ByVal};
				EigenType H = {HxVal, HyVal};
				sequenceVals.push_back(std::pair<EigenType, EigenType>(H,B));
			}
			MeasurementSequences.push_back(sequenceVals);

			sequence = sequence.NextSibling();

			sequenceNo +=1;



		}


	}

	void construct_measurement_matrices()
	{
		std::size_t numberOfSelectedPoints = 0;
		std::size_t numberOfSelectedSequences = 0;
		if (sequenceSelection.size() == 0)
		{
			numberOfSelectedPoints = numberOfDataPoints;
			numberOfSelectedSequences = MeasurementSequences.size();
		}
		else
		{
			for (std::size_t i = 0; i < sequenceSelection.size(); i++)
			{
				if (sequenceSelection[i] >  MeasurementSequences.size()-1)
				{
					std::cout << "ERROR: Sequence " << sequenceSelection[i] + 1  << " was selected " << " only " << MeasurementSequences.size() << " sequences given.";
					exit(1);
				}
				numberOfSelectedPoints += MeasurementSequences[sequenceSelection[i]].size();
			}
			numberOfSelectedSequences = sequenceSelection.size();
		}


		const int noOfHysterons = identModel.get_number_of_hysterons();
		A.resize(2*numberOfSelectedPoints, noOfHysterons);
		b.resize(2*numberOfSelectedPoints,1);
		std::size_t rowIndex = 0;


		for (std::size_t i = 0; i < numberOfSelectedSequences; i++)
		{
			int select = sequenceSelection[i];
			// TODO: Get rid of that


			EigenType dummyVec = MeasurementSequences[select][0].first *-1;
			identModel.calculate_evolution_H(dummyVec);
			identModel.accept_evolution();

			for (std::size_t j = 0; j < MeasurementSequences[select].size(); j++)
			{
				Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ABlock;
				ABlock.resize(2 ,noOfHysterons);
				identModel.gather_unit_magnetizations(MeasurementSequences[select][j].first, ABlock);
				EigenType J = MeasurementSequences[select][j].second - MeasurementSequences[select][j].first * constants::mu0;

				A.block(rowIndex,0,2, noOfHysterons) = ABlock;
				b.block(rowIndex,0,2,1) = J;
				rowIndex += 2;
			}

			identModel.reset_all_hysterons();
		}

	}

	void construct_regularization_matrices()
	{

	}

	void construct_matrices()
	{
		std::cout << "Building measurement matrices;\n";
		construct_measurement_matrices();
		std::cout << "Completed measurement matrices;\n";
		std::cout << "Resulting dimensions: ";
		std::cout << " A: "<< A.rows() << "x" << A.cols();
		std::cout << " b: "<< b.rows() << "x" << b.cols() <<std::endl;
	}

	void write_matrices_to_hdf5(const std::string& matrixPath)
	{
		const H5std_string filename (ewt_tools::concatenate_folders(matrixPath.c_str(),"identMatrices.h5"));
		std::cout << "Storing identification matrices in: " << filename << "\n";

		int rows;
		int cols;
		int rank;
		hsize_t     dimsf[2];


		H5::H5File file(filename, H5F_ACC_TRUNC);


		rows = b.rows();
		cols = b.cols();
		rank = 2;
		dimsf[0] = rows;
		dimsf[1] = cols;

		H5::DataSpace dataspaceb(rank, dimsf);
		H5::DataSet datasetb = file.createDataSet("b", H5::PredType::NATIVE_DOUBLE, dataspaceb);
		datasetb.write(b.data(), H5::PredType::NATIVE_DOUBLE);


		rows = A.rows();
		cols = A.cols(); // dataset dimensions
		rank = 2;
		dimsf[0] = rows;
		dimsf[1] = cols;
		H5::DataSpace dataspaceA(rank, dimsf);

		H5::DataSet datasetA = file.createDataSet("A", H5::PredType::NATIVE_DOUBLE, dataspaceA);
		datasetA.write(A.data(), H5::PredType::NATIVE_DOUBLE);


	}


};



#endif /* HYSTERESIS_OPERATOR_IDENTIFICATION_DPC_IDENTIFICATION_2D_H_ */
