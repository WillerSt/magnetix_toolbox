// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once
#include "magnetics_toolbox/mag_tools_basic.h"
#include "magnetics_toolbox/constants.h"
#include "magnetics_toolbox/quadrature_dof_coupler.h"
#include "magnetics_toolbox/scenario_description.h"
#include "magnetics_toolbox/mesh_container.h"
#include <Eigen/Core>

namespace mag_tools{
    class domain_rotation{
        public:
        using T = PetscScalar;
        using EigenVec = Eigen::Matrix<double, 2, 1>;
        using EigenMat = Eigen::Matrix<double, 2, 2>;
        const std::shared_ptr<dolfinxMesh<T>> mesh;
        const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> meshMarkers;
        const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> boundaryMarkers;

        const std::vector<size_t> rotIdx;
        const EigenVec rotCenter;
        const std::vector<std::vector<int>> compIndices;
        const std::vector<EigenVec> originalCoords;
        std::vector<EigenVec> currentCoords = originalCoords;


        const int masterInterfaceIdx;
        const int slaveInterfaceIdx;
  
        std::vector<double> rotationAngles;

        std::vector<mag_tools::quadrature_dof_coupler<T>> rotationCoupling;

        double lastAngleRad = 0.0;
        double lastAngleDeg = 0.0;

        public:
        domain_rotation(const mesh_container<T>& meshContainer,
             const EigenVec& rotCenterIn, 
             const std::vector<size_t>& rotIdxIn, 
             const int& masterIdx, const int& slaveIdx):
        mesh(meshContainer.get_mesh()), meshMarkers(meshContainer.get_mesh_markers()), 
        boundaryMarkers(meshContainer.get_boundary_markers()), rotIdx(rotIdxIn),rotCenter(rotCenterIn),
        compIndices(get_geometry_indices()), originalCoords(get_orig_coords()), 
        masterInterfaceIdx(masterIdx), slaveInterfaceIdx(slaveIdx) {

        }

        domain_rotation(const mesh_container<T>& meshContainer,  const mag_tools::scen::xml_model_entity& param):
            domain_rotation(meshContainer, {scen::as_double(param, "rotCenter_x"), scen::as_double(param, "rotCenter_y")}, 
            scen::determine_index_vector(param), scen::as_int(param, "masterInterface"),scen::as_int(param, "slaveInterface")){
                set_angle_vector(param);
        }
        void perform_rotation(const double& angle){
            double angRad = angle/180*constants::pi;
            EigenMat rotMat;
            rotMat << std::cos(angRad) , -std::sin(angRad) , std::sin(angRad) ,  std::cos(angRad); 
            

            for (std::size_t i = 0; i < currentCoords.size(); i++){
                currentCoords[i]= rotMat*(originalCoords[i]-rotCenter) + rotCenter;
                this->mesh->geometry().x()[compIndices[i][0]] = currentCoords[i][0];
                this->mesh->geometry().x()[compIndices[i][1]] = currentCoords[i][1];
            }

            lastAngleDeg = angle;
            lastAngleRad = angRad;         
        }
        void perform_rotation(const size_t& i){
            if (i>=rotationAngles.size()){
                std::cout << "ERROR: index " << i << " exceeds " << " size of defined angles " << rotationAngles.size() << std::endl;
                exit(1);  
            }
            perform_rotation(rotationAngles[i]);
        }
        void reverse_rotation(){
            for (std::size_t i = 0; i < originalCoords.size(); i++){
                this->mesh->geometry().x()[compIndices[i][0]] = originalCoords[i][0];
                this->mesh->geometry().x()[compIndices[i][1]] = originalCoords[i][1];
            }   
        }

        double get_last_angle_deg() const{
            return lastAngleDeg;
        }
        
        double get_last_angle_rad() const{
            return lastAngleRad;
        }
        void initialize_projection_coupling(const std::shared_ptr<dolfinxFunction<T>>& rotFunc){
            for (auto& idx: this->rotIdx){
                this->rotationCoupling.push_back(mag_tools::quadrature_dof_coupler<T>(rotFunc, this->meshMarkers, idx, this->mesh->topology()->dim(), this->mesh->topology()->dim()));
            }
            
        }

        void perpare_result_projection(const std::size_t& i){
            double angRad = this->rotationAngles[i]/180*constants::pi;
            EigenMat rotMat;
            rotMat << std::cos(angRad) , -std::sin(angRad) , std::sin(angRad) ,  std::cos(angRad); 
            for (auto& cplg: this->rotationCoupling){

                auto couplingMatrices = cplg.get_ref_vec();
                for (auto& fMat: *couplingMatrices){
                    *(fMat[0][0]) = rotMat.coeff(0,0);
                    *(fMat[1][0]) = rotMat.coeff(1,0);
                    *(fMat[1][1]) = rotMat.coeff(1,1);
                    *(fMat[0][1]) = rotMat.coeff(0,1);
                }
            }
        }
        
        void set_angle_vector(const mag_tools::scen::xml_model_entity& parameters){
            auto seriesXML = scen::as_string(parameters, "xmlFile");
                
            if (!boost::filesystem::exists(seriesXML)){
                std::cout << "ERROR: Time series XML " << seriesXML << " does not exist, aborting\n";
                exit(1);
            }
            else{
                rotationAngles = scen::read_xml_series(seriesXML, "rotationAngle");
                
                std::cout << "Read rotation angles:";
                for (std::size_t i = 0; i < rotationAngles.size(); i++){
                    std::cout <<" " << rotationAngles[i];
                }
                std::cout << std::endl;           
            }
        }

        std::vector<double> get_angle_vec(){
            return this->rotationAngles;
        }

        private:
        std::vector<std::vector<int>> get_geometry_indices(){

            auto coordFS = mesh->geometry().x();
            std::vector<int32_t> cells;
            for(auto& idx:rotIdx){
                auto partCells = meshMarkers->find(idx);
                cells.insert(cells.end(), partCells.begin(), partCells.end());
            }
            //mesh->topology()->create_connectivity(mesh->topology()->dim(), 0); // TODO: move to mesh_container
            std::vector<int> rotIndeces;
            for (auto& cell:cells){
                auto cellVertices = mesh->topology()->connectivity(mesh->topology()->dim(),0)->links(cell);
                for (auto& idx: cellVertices){
                    rotIndeces.push_back(idx);
                }
            }
            std::sort(rotIndeces.begin(), rotIndeces.end());
            rotIndeces.erase(std::unique(rotIndeces.begin(), rotIndeces.end()), rotIndeces.end());
            
            std::vector<std::vector<int>> coordIndices;
            for (auto& idx : rotIndeces){
                coordIndices.push_back({idx*3, idx*3+1, idx*3+2});                  
            } 
            return coordIndices;

        }

        std::vector<EigenVec> get_orig_coords(){
            std::vector<EigenVec> originalCoords;
            for (std::size_t i = 0; i<this->compIndices.size();i++){
                originalCoords.push_back({mesh->geometry().x()[compIndices[i][0]], mesh->geometry().x()[compIndices[i][1]]});
            }
            return originalCoords;                       
        }
    };

}