#include <magnetics_toolbox/scenario_description.h>
namespace mag_tools::scen{
    std::size_t search_for_name(const std::vector<mag_tools::scen::xml_parameter_node>& parameters,const std::string& name, const bool& abortIfMissing){
        for (std::size_t i = 0; i<parameters.size(); i ++){
            if (std::strcmp(name.c_str(), parameters[i].name.c_str() )== 0){
                return i;
            }
        }
        std::cout << "Error: Paramter with name " << name << "not found\n"; 
        for (std::size_t i = 0; i<parameters.size(); i ++){  
            parameters[i].print_parameter();
        }

        if (abortIfMissing == true){
            exit(1);
        }
        else{
            return std::size_t(-1);
        }         
    }

    std::vector<std::string> search_for_vector(const std::vector<mag_tools::scen::xml_parameter_node>& parameters, const std::string& name, const std::size_t& startIndex){
        std::vector<std::string> stringVec;
        std::size_t idx = startIndex;
        for (std::size_t i = 0; i<parameters.size(); i ++){
            std::string s = name;
            s += "_";
            s += std::to_string(idx);
            if (std::strcmp(s.c_str(), parameters[i].name.c_str() )== 0){
                stringVec.push_back(s);
                idx += 1;
            }
            
        }
        return stringVec;
    }

    std::vector<double> read_xml_series(const std::string& seriesXML, const std::string& quantity){
        std::vector<double> resultVec;
        pugi::xml_document doc;
        if (!boost::filesystem::exists(seriesXML)){
            std::cout << "ERROR: Time series XML" << seriesXML << " does not exist, aborting\n";
            exit(1);
        }
        else{
            doc.load_file(seriesXML.c_str());
            std::stringstream reader;
            auto timeSeries = doc.child(DEFINE_TIME_SERIES);
            for (auto node = timeSeries.child(DEFINE_TIME_VARIATION); node; node=node.next_sibling(DEFINE_TIME_VARIATION)){
                if (strcmp(node.attribute("quantity").as_string(), quantity.c_str())==0){
                    reader << node.text().as_string();
                    break;
                }
            }
            double value;
            while (reader >> value){
                resultVec.push_back(value);
            }
        }
        return resultVec;
    }

    std::vector<double> read_xml_material(const std::string& seriesXML, const std::string& quantity){
        std::vector<double> resultVec;
        pugi::xml_document doc;
        if (!boost::filesystem::exists(seriesXML)){
            std::cout << "ERROR: Material XML" << seriesXML << " does not exist, aborting\n";
            exit(1);
        }
        else{
            doc.load_file(seriesXML.c_str());
            std::stringstream reader;
            auto timeSeries = doc.child(DEFINE_MATERIAL_DEFINITION);
            for (auto node = timeSeries.child(DEFINE_MATERIAL_DATA); node; node=node.next_sibling(DEFINE_MATERIAL_DATA)){
                if (strcmp(node.attribute("quantity").as_string(), quantity.c_str())==0){
                    reader << node.text().as_string();
                    break;
                }
            }
            double value;
            while (reader >> value){
                resultVec.push_back(value);
            }
        }
        return resultVec;
    }

    std::string xml_model_entity::get_type_parameter() const{
            std::size_t idx = search_for_name(this->parameter, "type", false);
            if (idx == std::size_t(-1)){
                return "";
            }
            else{
                return this->parameter[idx].value;
            }
            
        }


    double as_double(const std::vector<mag_tools::scen::xml_parameter_node>& parameters, const std::string& name){
        return  std::stod(parameters[scen::search_for_name(parameters, name)].value.c_str());
    }
    double as_double(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name){
        auto parameters = entityDesc.get_parameters();
        return as_double(parameters, name);
    }



    int as_int(const std::vector<mag_tools::scen::xml_parameter_node>& parameters, const std::string& name){
        return  std::stoi(parameters[scen::search_for_name(parameters, name)].value.c_str());
    }

    int as_int(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name){
        auto parameters = entityDesc.get_parameters();
        return as_int(parameters, name);
    }

    std::vector<std::size_t> determine_index_vector(const mag_tools::scen::xml_model_entity& params){
        auto domainVec = scen::search_for_vector(params.get_parameters(), "domain");
        std::vector<std::size_t> idxVec;
        for (auto& dom: domainVec){
            idxVec.push_back(scen::as_int(params, dom));
        }
        return idxVec;
    }

    

    std::string as_string(const std::vector<mag_tools::scen::xml_parameter_node>& parameters, const std::string& name){
        return  parameters[scen::search_for_name(parameters, name)].value;
    }

    std::string as_string(const mag_tools::scen::xml_model_entity& entityDesc, const std::string& name){
        auto parameters = entityDesc.get_parameters();
        return as_string(parameters, name);
    }

}