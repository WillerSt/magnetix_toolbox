// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once

#include <pugixml.hpp>
#include <string>

#include <iostream>
#include <boost/filesystem.hpp>
#include <pugixml.hpp>

#include <petscsys.h>


namespace mag_tools{
    constexpr auto DEFINE_MESH_DEFINITION = "MeshDefinition";
    constexpr auto DEFINE_SOLVER_PARAMETERS = "SolverParameters";
    constexpr auto DEFINE_BOUNDARY_DEFINTIION = "BoundaryhDefinition";
    constexpr auto DEFINE_MATERIAL_DEFINITION = "MaterialDefinition";
    constexpr auto DEFINE_SOURCE_DEFINITION = "SourceDefinition";
    constexpr auto DEFINE_ENTITY_DEFINITION = "EntityDefintion";
    constexpr auto DEFINE_SCENARIO_DESCRIPTION = "ScenarioDescription";
    constexpr auto DEFINE_TIME_STEPPING = "TimeStepping";
    constexpr auto DEFINE_TIME_VARIATION = "TimeVariation";
    constexpr auto DEFINE_TIME_SERIES = "TimeSeries";
    constexpr auto DEFINE_FORCE_CALCULATION = "ForceCalculation";
    constexpr auto DEFINE_MOVING_DOMAIN = "MovingDomain";
    constexpr auto DEFINE_MATERIAL_DATA = "MaterialData";
    constexpr auto DEFINE_POINT_EVAL = "PointEvaluation";

    namespace scen{
        class mesh_input{
            public:
            const std::string meshFile;
            const std::string bdryFile;

            public:
            mesh_input(const std::string& meshFileIn, const std::string& bdryFileIn):
            meshFile(meshFileIn), bdryFile(bdryFileIn){
                if (!boost::filesystem::exists(meshFile)){
                    std::cout << "ERROR: Mesh file " << meshFile << "does not exist" << std::endl;
                    exit(1);
                }
                if (!boost::filesystem::exists(bdryFile)){
                    std::cout << "ERROR: Boundary file " << bdryFile << "does not exist" << std::endl;
                    exit(1);
                }
            }

            std::string mesh_file() const{
                return this->meshFile;
            }
            std::string boundary_file() const{
                return this->bdryFile;
            }


        };

        class subdomain_description{
            public:
                const std::string materialType;
                const std::size_t sudomainIndex;
        };

        class xml_parameter_node{
        public:
        const std::string name;
        std::string value;
        const std::string description;

        xml_parameter_node(const std::string& nameIn, const std::string& valueIn, const std::string& descriptionIn = ""):
            name(nameIn), value(valueIn), description(descriptionIn){}

        public: 
        void print_parameter() const{
            std::cout << "name: \"" << name << "\";" << " value: \"" << value << "\";"; 
            if (!description.empty()){
                std::cout << " desc: \"" << description << "\"";
            }
            std::cout<< std::endl;
        }
    };

    class xml_model_entity{

        public:
        const std::string type;
        
        private:
        std::vector<xml_parameter_node> parameter;

        public:
        xml_model_entity(const std::string& typeIn): type(typeIn){

        }
        void add_parameter(const xml_parameter_node& singleParam){
            parameter.push_back(singleParam);
        }

        std::vector<xml_parameter_node> get_parameters() const{
            return parameter;
        }
        void print_entity_info() const{
            std::cout << "\nEntity Type: \"" << type << "\""<< std::endl;
            for(auto& param: parameter){
                param.print_parameter();
            } 
        }
        std::string get_type_parameter() const;
        void check_files(const boost::filesystem::path& scenPath){

            int rank;
            MPI_Comm_rank( MPI_COMM_WORLD, &rank);

            for (auto& param: this->parameter){
                if ((param.name.find("file") != std::string::npos)||(param.name.find("File") != std::string::npos)){
                    std::string searchFile = param.value.c_str();
                    if (!boost::filesystem::exists(searchFile)){
                        auto candidateFile1 = scenPath / searchFile;
                        auto candidateFile2 = scenPath / "input"/ searchFile;
                        if (boost::filesystem::exists(candidateFile1)){
                            searchFile = candidateFile1.c_str();
                        }
                        else if (boost::filesystem::exists(candidateFile2)){
                            searchFile = candidateFile2.c_str();
                        }
                        else{
                            if (rank == 0){
                                std::cout << "ERROR: Could not find file " << param.name <<": " << searchFile <<std::endl;
                            }
                           exit(1);
                        }
                        param.value = searchFile;
                    }
                    if (rank == 0){
                        std::cout << "Found " << param.name <<": " << searchFile << std::endl;
                    }
                }
            }

        }
    };

    std::size_t search_for_name(const std::vector<mag_tools::scen::xml_parameter_node>& parameters,const std::string& name, const bool& abortIfMissing = true);

    std::vector<std::string> search_for_vector(const std::vector<mag_tools::scen::xml_parameter_node>& parameters, const std::string& name, const std::size_t& startIndex = 1);

    std::vector<double> read_xml_series(const std::string& seriesXML, const std::string& quantity);
    std::vector<double> read_xml_material(const std::string& seriesXML, const std::string& quantity);
    
    double as_double(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name);
    int as_int(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name);
    std::string as_string(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name);
    std::vector<std::size_t> determine_index_vector(const mag_tools::scen::xml_model_entity& params);

    }





    class scenario_description{
        public:
        const std::string xmlFilename;

        private:
        const pugi::xml_document doc;
        const boost::filesystem::path scenPath;

        using entityVec = std::vector<scen::xml_model_entity>;

        public:
        const std::string scenarioName;
        const mag_tools::scen::mesh_input meshInput;
        const std::vector<scen::xml_model_entity> materialInput;
        const std::vector<scen::xml_model_entity> sourceInput;
        const scen::xml_model_entity solverInput;
        const scen::xml_model_entity timeStepping;
        const std::vector<scen::xml_model_entity> boundaryDef;
        const entityVec forceCalc;
        const entityVec movingDomains;
        const entityVec pointEval;
        const std::vector<std::string> outputConfig;



        public:
        scenario_description(const std::string& xmlFilenameIn): xmlFilename(xmlFilenameIn), doc(load_document()), 
            scenPath(boost::filesystem::canonical(boost::filesystem::path(this->xmlFilename).parent_path())),
            scenarioName(get_scenario_name()),
            meshInput(get_mesh_input()), 
            materialInput(get_entitiy_group(DEFINE_MATERIAL_DEFINITION)), 
            sourceInput(get_entitiy_group(DEFINE_SOURCE_DEFINITION)), 
            solverInput(get_entitiy_group(DEFINE_SOLVER_PARAMETERS)[0]),
            timeStepping(get_entitiy_group(DEFINE_TIME_STEPPING)[0]),
            boundaryDef(get_entitiy_group(DEFINE_BOUNDARY_DEFINTIION)),
            forceCalc(get_entitiy_group(DEFINE_FORCE_CALCULATION)),
            movingDomains(get_entitiy_group(DEFINE_MOVING_DOMAIN)),
            pointEval(get_entitiy_group(DEFINE_POINT_EVAL)),
            outputConfig({"not implemented yet"}){

        }

        scen::xml_model_entity get_material_model_input(const std::size_t& index) const{
            return materialInput[index];
        }
        std::string scenario_name() const{
            return scenarioName;
        }

        std::string get_scen_path(){
            return this->scenPath.c_str();
        }

        private:


            std::vector<scen::xml_model_entity> get_entitiy_group(const std::string& groupName){
                
                std::vector<scen::xml_model_entity> entityVec;
                for (pugi::xml_node tool = doc.child(DEFINE_SCENARIO_DESCRIPTION).child(groupName.c_str()); tool; tool = tool.next_sibling(groupName.c_str())){     
                    scen::xml_model_entity groupMember(groupName);  

                    for (pugi::xml_node param = tool.child(DEFINE_ENTITY_DEFINITION); param; param = param.next_sibling(DEFINE_ENTITY_DEFINITION)){
                        groupMember.add_parameter(scen::xml_parameter_node(param.attribute("name").value(), param.attribute("value").value()));
                    }
                    groupMember.check_files(this->scenPath);
                    entityVec.push_back(groupMember);

                }
                return entityVec;

            }

            pugi::xml_document load_document(){
                pugi::xml_document docProto;

                if (!boost::filesystem::exists(xmlFilename)){
                    std::cout << "Scenario file " << xmlFilename << " does not exist, aborting execution" << std::endl;
                    exit(1);
                }
                else{
                    docProto.load_file(this->xmlFilename.c_str());
                }
                return docProto;
            }

            std::string get_scenario_name(){
                std::string name = doc.child(DEFINE_SCENARIO_DESCRIPTION).attribute("name").as_string();
                if (name.empty()){
                    std::cout << "Name attribute in " << DEFINE_SCENARIO_DESCRIPTION << " required \n";
                    exit(1);
                }

                int rank;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);

                if (rank == 0){
                    std::cout << "Scenario file: " << xmlFilename << std::endl;
                    std::cout << "New top folder:" << scenPath.c_str() << std::endl;
                } 

                return doc.child(DEFINE_SCENARIO_DESCRIPTION).attribute("name").as_string();

            }


            std::vector<scen::xml_model_entity> get_source_description(){
                std::vector<scen::xml_model_entity> srcVec;
                return srcVec;

            }
            scen::mesh_input get_mesh_input(){
                auto meshEntities = get_entitiy_group(DEFINE_MESH_DEFINITION);
                if (meshEntities.size() > 1){
                    std::cout << "Warning: There should not be more than one " << DEFINE_MESH_DEFINITION << " node, defaulting to first entry" << std::endl;
                }
                auto idxMesh = scen::search_for_name(meshEntities[0].get_parameters(), "meshFile");
                auto idxBdry= scen::search_for_name(meshEntities[0].get_parameters(), "boundaryFile");

                auto meshFile = meshEntities[0].get_parameters()[idxMesh].value;
                auto bdryFile = meshEntities[0].get_parameters()[idxBdry].value;
                return mag_tools::scen::mesh_input(meshFile, bdryFile);
                
                
            }
            


    };
}