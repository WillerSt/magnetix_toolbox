#include <magnetics_toolbox/file_management.h>
#include <fstream>
#include <ctime>
#include <iostream>
#include <dolfinx.h>
#include <boost/program_options.hpp>

namespace mag_tools{


	std::string get_scenario_file(int argc, char* argv[]){
		boost::program_options::options_description desc{"Options"};
		desc.add_options()
			("scen, Scen, scenario, Scenario", boost::program_options::value<std::string>()->default_value(""), "Scenario file");

		boost::program_options::variables_map inputParameter;
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), inputParameter);
		return inputParameter["scen"].as<std::string>();
	}

	file_manager::file_manager(const std::string& topFolder, const std::string& scenarioName, const bool& addDate){
		int rank;
		MPI_Comm_rank( MPI_COMM_WORLD, &rank);
			scenarioPath = boost::filesystem::path (topFolder);

		if (addDate)
		{
			char date[9];
			std::time_t t=std::time(NULL);
			std::strftime(date, sizeof(date), "%Y%m%d", std::localtime(&t));
			scenarioPath /= scenarioName+"_"+date;
		}
		else
		{
			scenarioPath /= scenarioName;
		}

		resultPath = scenarioPath / "results";
		inputPath = scenarioPath / "input"; 
		logPath = scenarioPath / "logs";

		if (rank == 0){
			this->create_folder_structure(topFolder);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}

    void file_manager::create_folder_structure(const std::string& topFolder)
    {
        if (!boost::filesystem::exists(topFolder))
        {
        boost::filesystem::create_directory(topFolder);
        }
		
		auto parPath = scenarioPath.parent_path();

        boost::filesystem::create_directory(scenarioPath);
        boost::filesystem::create_directory(inputPath);
        boost::filesystem::create_directory(resultPath);
        boost::filesystem::create_directory(logPath);

    }

    void file_manager::gather_file_entries(const std::vector<mag_tools::scen::xml_model_entity>& entityVec, std::vector<std::string>& filesToCopy) const{
        for (auto& srcIn: entityVec){
            for(auto& param: srcIn.get_parameters()){
                if ((param.name.find("file") != std::string::npos)||(param.name.find("File") != std::string::npos)){
                    filesToCopy.push_back(param.value);
                    //std::cout << "found file " << filesToCopy.back() << std::endl;
                    continue;
                };
            }
        }
    }

    void file_manager::copy_all_input_files(const mag_tools::scenario_description& scenIn) const{
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        
        if (rank == 0){
            std::vector<std::string> filesToCopy;
            //filesToCopy.push_back(scenIn.xmlFilename);
            //filesToCopy.push_back(scenIn.meshInput.boundary_file());
            //filesToCopy.push_back(scenIn.meshInput.mesh_file());
            gather_file_entries(scenIn.sourceInput, filesToCopy);
            gather_file_entries(scenIn.materialInput, filesToCopy);
           this->copy_input_files(filesToCopy);
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }

  

  void file_manager::copy_input_files(const std::vector<std::string>& inputFiles) const{

    for (std::size_t i = 0; i<inputFiles.size(); i++)
    {
      boost::filesystem::path filePath(inputFiles[i]);
      boost::filesystem::path destFile = this->inputPath /filePath.filename() ;
      // boost::filesystem::path sourceFile = boost::filesystem::path(".") / inputFiles[i];
      boost::filesystem::path sourceFile = boost::filesystem::path( inputFiles[i]);

      boost::filesystem::path p(destFile);
      boost::filesystem::create_directories(p.parent_path());
      //std::cout << "Copying file " << sourceFile.c_str() << " to " << destFile.c_str() << std::endl;
      boost::filesystem::copy_file(sourceFile, destFile, boost::filesystem::copy_option::overwrite_if_exists);
    }
  }

  void file_manager::write_matlab_import_file()
  {
    std::string mFile ((this->scenarioPath/"importResults.m").c_str());
    std::ofstream outStream;
    outStream.open(mFile);
    outStream << "pathToResults = sprintf(\'%s/results\', pwd); \n"
              << "mkdir('matlab'); \n"
              << "results = vtk_data_reader(pathToResults,[],sprintf(\'%s/matlab\',pwd)); \n";
    outStream.close();
  }
}