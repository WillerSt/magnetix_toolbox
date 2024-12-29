/*
* Creates a folder structure equal to:
*
*                                 |--results
* topFolder--scenarioName(_date)--|--input
*                                 |--logs
*
* name of subfolders is mirrored in vtk-output etc.
*/
#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <magnetics_toolbox/scenario_description.h>

namespace mag_tools{
    std::string get_scenario_file(int argc, char* argv[]);
    

    class file_manager{

        private:
         boost::filesystem::path scenarioPath;
         boost::filesystem::path resultPath;
         boost::filesystem::path inputPath;
         boost::filesystem::path logPath;
         

        public:
        file_manager(const std::string& topFolder, const std::string& scenarioName, const bool& addDate);

        void copy_all_input_files (const mag_tools::scenario_description& scenIn) const;

        const boost::filesystem::path scenario_path() const {
            return scenarioPath;
        }
        const boost::filesystem::path result_path() const {
            return resultPath;
        }
        const boost::filesystem::path input_path() const {
            return inputPath;
        }
        const boost::filesystem::path log_path() const {
            return logPath;
        }

        public:
        static void abort_if_not_existent (const std::string& file, const std::string& requFrom){
            if (!boost::filesystem::exists(file)){
                std::cout << "File required from " << requFrom << " does not exist:\n" << file << std::endl;
                exit(1); 
            }
        }
        
        private:
        void create_folder_structure(const std::string& topFolder);
        void copy_input_files(const std::vector<std::string>& inputFiles) const;
        void write_matlab_import_file();
        void gather_file_entries(const std::vector<mag_tools::scen::xml_model_entity>& entityIn, std::vector<std::string>& filesToCopy) const;

    };
}

