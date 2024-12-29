#pragma once
#include <vector>
#include <pugixml.hpp>
#include <dolfinx.h>
#include <iostream>


namespace mag_tools{
    namespace output{
        class time_series{
            public:
            using EigenVec = Eigen::Matrix<double, 2, 1>;
            
            const std::string quantity;
            std::vector<double> values;
            std::vector<std::string> attrNames;
            std::vector<std::string> attrValues;

            public:
            time_series(const std::string& quantityName, const std::vector<std::string>& attributesIn={}, const std::vector<std::string>&attrValusIn = {}):
                quantity(quantityName), attrNames(attributesIn), attrValues(attrValusIn){

            }
            void add_value(const double& newRes){
                this->values.push_back(newRes);
            }

            void add_values(const std::vector<double>& newRes){
                for (size_t i = 0; i < newRes.size(); i++){
                    this->values.push_back(newRes[i]);
                }                 
            }
            void add_values(const std::vector<EigenVec> newRes){
                for (size_t i = 0; i < newRes.size(); i++){
                    this->values.push_back(newRes[i][0]);
                    this->values.push_back(newRes[i][1]);
                }  
            }

        };

        class result_xml{
            std::vector<std::shared_ptr<time_series>> results;

            pugi::xml_document xmlDoc = pugi::xml_document();
            pugi::xml_node desc = xmlDoc.append_child(pugi::node_declaration);
            pugi::xml_node root = xmlDoc.append_child("TimeSeries");

            const std::string resultFile;

            const int rank;

            public:
            result_xml(const std::string& resultFileIn):resultFile(resultFileIn), rank(get_rank()){
                desc.append_attribute("version") = "1.0";
                
            }
            void write_file(){
                for (auto& series: results){
                    
                    auto seriesNode = root.append_child("TimeVariation");
                    seriesNode.append_attribute("name") = series->quantity.c_str();
                    for (size_t i = 0; i<series->attrNames.size(); i++){
                        seriesNode.append_attribute(series->attrNames[i].c_str()) = series->attrValues[i].c_str();
                    }
                    auto textContent = seriesNode.text();
                    std::stringstream textString;
                    for (size_t i  = 0; i < series->values.size(); i++){
                        textString << " " << series->values[i];
                    }
                    textContent.set(textString.str().c_str());
                }


                if (rank == 0){
                    xmlDoc.save_file(resultFile.c_str());
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }

            std::shared_ptr<time_series> append_time_series(const std::string& name, const std::vector<std::string>& attributesIn={}, const std::vector<std::string>& attrValusIn = {}){
                auto timeSeries = std::make_shared<time_series>(name, attributesIn, attrValusIn);
                results.push_back(timeSeries);
                return timeSeries;
            }

            private:
            int get_rank(){
                int myRank;
                MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
                return myRank;
            }
        };
    }
}
