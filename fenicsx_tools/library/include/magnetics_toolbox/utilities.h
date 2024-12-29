#pragma once
#include "pugixml.hpp"
#include <vector>
#include <petsc.h>

namespace mag_tools{
    namespace utils{
        class iteration_monitor{
            std::vector<int> seriesSteps = {};
            std::vector<std::vector<double>> iterations = {};


            pugi::xml_document xmlDoc = pugi::xml_document();
            pugi::xml_node desc = xmlDoc.append_child(pugi::node_declaration);
            pugi::xml_node iterStatNode = xmlDoc.append_child("IterationStatistic");

            const int rank;

            public:
            iteration_monitor():rank(get_rank()){
                desc.append_attribute("version") = "1.0";
            }
            void next_step(){
                iterations.push_back({});
            }
            void next_iteration(const double& res){
                iterations.back().push_back(res);
            }
            void write_statistic_file(const std::string& fileName){
                 if (rank == 0){
                    for (size_t step = 0; step < iterations.size(); step++){
                        auto iterationsNode = iterStatNode.append_child("Iteration");
                        iterationsNode.append_attribute("step") = step;
                        auto textContent = iterationsNode.text();
                        std::stringstream textString;
                        for (size_t i  = 0; i < iterations[step].size(); i++){
                            textString << " " << iterations[step][i];
                        }
                        textContent.set(textString.str().c_str());
                    }                               
                    xmlDoc.save_file(fileName.c_str());
                }
                MPI_Barrier(MPI_COMM_WORLD);              
            }

            private:
            int get_rank() const{
                int currentRank;
                MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);
                return currentRank;
            }
        };

    }    
}