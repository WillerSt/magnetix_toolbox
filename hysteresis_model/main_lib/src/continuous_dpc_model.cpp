// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include "Hysteresis_Operator/continuous_dpc_model.h"

#ifdef ENALBE_DPC_CONT_TIMING
Timer dpcContTimer;
#endif


template<int dimension>
std::vector<std::shared_ptr<const CS_sphere_cont<dimension>>> continuous_dpc_model<dimension>::initialize_cont_spheres()
{
	std::vector<std::shared_ptr<const CS_sphere_cont<dimension>>> result;
	for (std::size_t i = 0; i < (*mallocGrid).numberOfSurfaces; i++)
	{
			result.push_back(this->mallocGrid->get_pointer_to_surface(i));
			this->numberOfHysterons += 1;

			if (i != (*result[i]).get_surfaceID())
			{
				std::cout << "WARNING: Indices of hysterons are not equal to IDs of critical surfaces \n";
			}
	}

	return result;
}

template<int dimension>
void continuous_dpc_model<dimension>::write_path_visualization_file(const std::string& fileNameRaw, const std::string& folderName) const
	{

		//std::string fileName(ewt_tools::concatenate_folders(folderName.c_str(),fileNameRaw.c_str()));
		std::string fileName = (boost::filesystem::path(folderName)/fileNameRaw).string(); 
		if (!boost::filesystem::exists(folderName.c_str()))
		{
			boost::filesystem::create_directories(folderName);
		}
		std::stringstream dataString1, dataString2, dataString3;
		tinyxml2::XMLDocument dataXML;
		auto decl = dataXML.NewDeclaration();
		dataXML.InsertEndChild(decl);

		auto rootNode = dataXML.NewElement("DPC_visual");

		auto evolNode = dataXML.NewElement("MagFieldEvolution");

		evolNode->SetAttribute("dim", 2);
		evolNode->SetAttribute("steps", int(Hevol.size()-1));

		for (std::size_t i = 0; i < Hevol.size(); ++i)
		{
			dataString1 << Hevol[i][0] << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< " ";
			dataString2 << Hevol[i][1] << std::setprecision(std::numeric_limits<double>::digits10 + 1)
			<< " ";
		}
		auto HxNode = dataXML.NewElement("Hx");
		auto HyNode = dataXML.NewElement("Hy");

		HxNode->SetText(dataString1.str().c_str());
		HyNode->SetText(dataString2.str().c_str());

		evolNode->InsertEndChild(HxNode);
		evolNode->InsertEndChild(HyNode);

		rootNode->InsertEndChild(evolNode);



		auto pathNode = dataXML.NewElement("ReducedPaths");
		pathNode->SetAttribute("points", int(combinedSurfaces->size()));

		for (std::size_t i = 0; i < (*combinedSurfaces).size(); ++i)
		{


			auto pointNode = dataXML.NewElement("Point");
			pointNode->SetAttribute("Index", int(i));

			dataString1.str(std::string());
			dataString1.clear();
			dataString1 << (*combinedSurfaces)[i].position.transpose() << std::setprecision(std::numeric_limits<double>::digits10 + 1);
			pointNode->SetAttribute("position", dataString1.str().c_str());

			dataString1.str(std::string());
			dataString1.clear();
			dataString1 << (*combinedSurfaces)[i].radMax << std::setprecision(std::numeric_limits<double>::digits10 + 1);
			pointNode->SetAttribute("radMax", dataString1.str().c_str());

			dataString1.str(std::string());
			dataString1.clear();
			dataString2.str(std::string());
			dataString2.clear();
			dataString3.str(std::string());
			dataString3.clear();

			EigenStdVector pathVector;
			build_relevant_path(pathVector, i);

			for (std::size_t j = 0; j < pathVector.size(); ++j)
			{
				dataString1 << pathVector[j][0] << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< " ";
				dataString2 << pathVector[j][1] << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< " ";
				dataString3 << (*combinedSurfaces)[i].distance(pathVector[j]) << std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< " ";
			}

			auto PxNode = dataXML.NewElement("Hx");
			auto PyNode = dataXML.NewElement("Hy");
			auto DistNode = dataXML.NewElement("distance");

			PxNode->SetText(dataString1.str().c_str());
			PyNode->SetText(dataString2.str().c_str());
			DistNode->SetText(dataString3.str().c_str());

			pointNode->InsertEndChild(PxNode);
			pointNode->InsertEndChild(PyNode);
			pointNode->InsertEndChild(DistNode);

			pathNode->InsertEndChild(pointNode);

		}

		rootNode->InsertEndChild(pathNode);
		dataXML.InsertEndChild(rootNode);
		dataXML.SaveFile(fileName.c_str());

	}


template class continuous_dpc_model<2>;
