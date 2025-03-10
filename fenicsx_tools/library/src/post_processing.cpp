// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#include "magnetics_toolbox/post_processing.h"
#include <magForce_virtualWork.h>
#include "magnetics_toolbox/utility.h"

namespace mag_tools{

    template<typename T> nodal_force<T>::nodal_force(const int& bdryIdx, const std::shared_ptr<dolfinxRHS<T>>& fVecIn,
                const std::shared_ptr<std::vector<double>>& dofCoords,
                const std::shared_ptr<const dolfinxFS<T>>& P1, 
                const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>>& boundaryMarkers,
                const EigenVec& rotCenterIn):
                idx(bdryIdx),
                fDofs(dolfinx::fem::locate_dofs_topological(*(P1->mesh()->topology_mutable()),
                    *(P1->dofmap()), P1->mesh()->topology()->dim()-1,boundaryMarkers->find(idx))),
                fVec(fVecIn), rotCenter(rotCenterIn){
            for (std::size_t i=0; i<fDofs.size(); i++){
                nodeCoords.push_back({(*dofCoords)[3*fDofs[i]], (*dofCoords)[3*fDofs[i]+1]});
                nodalForces.push_back({0.0, 0.0});
            }
        }


       template<typename T> void nodal_force<T>::update_force(const double& scaleFac, const double& angle){
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            double Fx = 0.0, Fy = 0.0, TzTemp = 0.0;

            for (std::size_t i = 0; i<fDofs.size(); i++){
                std::size_t idx = fDofs[i];
                
                nodalForces[i] = {fVec->array()[idx*2], fVec->array()[idx*2+1]};
                Fx = Fx + nodalForces[i][0];
                Fy = Fy + nodalForces[i][1];
                EigenVec lever = nodeCoords[i]-rotCenter;
                TzTemp  = TzTemp +  nodalForces[i][1]*lever[0] - nodalForces[i][0]*lever[1];
            }
            TzTemp =-TzTemp*scaleFac;
            Fy =-Fy*scaleFac;
            Fx =-Fx*scaleFac;
            
            

            // sum up actual force across all processes
            MPI_Allreduce(MPI_IN_PLACE, &Fx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Fy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &TzTemp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


            // update class variables
            F = {Fx, Fy};
            Tz = TzTemp;

            if (std::fabs(angle) > 1e-10){
                Eigen::Matrix<double, 2, 2> rotMat;
                rotMat << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);
                F = rotMat*F;
            }
            

            #ifdef __DEBUG__FMT__FORCE_OUTPUT
            if (rank == 0){ 
                    std::cout << "Fx " << F[0] <<std::endl;
                    std::cout << "Fy " << F[1] <<std::endl;
                    std::cout << "T " << Tz <<std::endl;

            }
            #endif //__DEBUG__FMT__FORCE_OUTPUT
            
        }

        template<typename T> void nodal_force<T>::store_results(const double& time){
            timeSeries.push_back(time);
            forceSeries.push_back(this->F);
            torqueSeries.push_back(this->Tz);
        }



    template<typename T>energy_density<T>::energy_density(const mesh_container<T>& meshContainer,
                    const std::shared_ptr<dolfinxFunction<T>>& BIn,
                    const std::shared_ptr<dolfinxFunction<T>>& HIn,
                    const std::shared_ptr<std::vector<mag_tools::material_container<T>>>matDefsIn ):
                   matDefs(matDefsIn), 
                   mesh(meshContainer.get_mesh()),
                   meshMarkers(meshContainer.get_mesh_markers()),
                   Wmag(std::make_shared<dolfinxFunction<T>>(mag_tools::utility::create_quadrature_functionspace<T>(mesh, "sca"))),
                   B(BIn), H(HIn){

                    Wmag->name = "Wmag0";
                    Wmag->x()->set(0.0);
                    for (auto& mat: *matDefs){
                        dofCouplers.push_back(std::make_shared<quadrature_dof_coupler<T>>(Wmag, meshMarkers, mat.get_domain_index(), 1, 1));
                        mat.couple_energy_density(dofCouplers[dofCouplers.size()-1]);
                    }
    }
    template<typename T> void energy_density<T>::perform_update(){
        for (size_t i = 0; i < WmagDofs.size(); i++){
            *WmagDofs[i][0][0] =  0.5*(*BDofs[i][0][0] * *HDofs[i][0][0] + *BDofs[i][1][0] * * HDofs[i][1][0]);
        }
        for (auto& mat: *matDefs){
            mat.update_energy_density();
        }
        
    }

    template<typename T> force_calculation<T>::force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                        const std::shared_ptr<dolfinxFunction<T>>& BIn,
                        const std::shared_ptr<dolfinxFunction<T>>& HIn,
                        const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                        const std::shared_ptr<energy_density<T>>& energyDensity,
                        const mag_tools::mesh_container<T>& meshContainer,
                        const std::vector<std::size_t>& idxIn,
                        const double& domDepthIn,
                        Eigen::Matrix<double, 2, 1> rotCenter):
                        A(AIn), B(BIn), H(HIn), nuDiff(nuDiffIn),
                        mesh(meshContainer.get_mesh()), boundaryMarkers(meshContainer.get_boundary_markers()), 
                        P1(std::make_shared<const dolfinxFS<T>>(
                            dolfinx::fem::create_functionspace(
                                mesh,
                                basix::create_element<U<T>>(
                                    basix::element::family::P, 
                                    dolfinx::mesh::cell_type_to_basix_type(mesh->geometry().cmap().cell_shape()), 
                                    1,
                                    basix::element::lagrange_variant::unset,
                                    basix::element::dpc_variant::unset, false),
                                {mesh->topology()->dim()})
                                )
                            ),
                        Wmag0(energyDensity->get_density_function()),
                        LForce(std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
                            *form_magForce_virtualWork_L, {P1}, {{"A", A}, {"B", B}, {"H", H}, {"nu", nuDiff}, {"Wmag0", Wmag0}, {"vacInd", vacInd}}, 
                                {   {"c0",std::make_shared<const fem::Constant<T>>(0.0)}, 
                                    {"c1",std::make_shared<const fem::Constant<T>>(0.0)}, 
                                    {"c2",std::make_shared<const fem::Constant<T>>(0.0)}}, {}, {}))),
                        domDepth(domDepthIn){

        for (auto& i : idxIn){
            forceCalculations.push_back(nodal_force<T>(i,this->F_vec, this->dofCoords, this->P1, this->boundaryMarkers, rotCenter));
        }
        
    }
    template<typename T> force_calculation<T>::force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                        const std::shared_ptr<dolfinxFunction<T>>& BIn,
                        const std::shared_ptr<dolfinxFunction<T>>& HIn,
                        const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                        const std::shared_ptr<energy_density<T>>& energyDensity,
                        const mag_tools::mesh_container<T>& meshContainer,
                        const mag_tools::scen::xml_model_entity& fCalcParam):
                        force_calculation(AIn, BIn, HIn, nuDiffIn,energyDensity, meshContainer,
                        determine_index_vector(fCalcParam), scen::as_double(fCalcParam, "domainDepth"), determine_rot_center(fCalcParam))
        {}

    template<typename T> force_calculation<T>::force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                    const std::shared_ptr<dolfinxFunction<T>>& BIn,
                    const std::shared_ptr<dolfinxFunction<T>>& HIn,
                    const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                    const std::shared_ptr<energy_density<T>>& energyDensity,
                    const mag_tools::mesh_container<T>& meshContainer,
                    const double& domDepthIn):A(AIn), B(BIn), H(HIn), nuDiff(nuDiffIn),
                        mesh(meshContainer.get_mesh()), meshMarkers(meshContainer.get_mesh_markers()), 
                        boundaryMarkers(meshContainer.get_boundary_markers()), 
                        P1(std::make_shared<const dolfinxFS<T>>(
                            dolfinx::fem::create_functionspace(
                                mesh,
                                basix::create_element<U<T>>(
                                    basix::element::family::P, basix::cell::type::triangle, 1,
                                    basix::element::lagrange_variant::unset,
                                    basix::element::dpc_variant::unset, false),
                                {2})
                                )
                            ),
                        Wmag0(energyDensity->get_density_function()),
                        vacInd(std::make_shared<dolfinxFunction<T>>(std::make_shared<const dolfinxFS<T>>(
                            dolfinx::fem::create_functionspace(
                                mesh,
                                basix::create_element<U<T>>(
                                    basix::element::family::P, basix::cell::type::triangle, 0,
                                    basix::element::lagrange_variant::unset,
                                    basix::element::dpc_variant::unset, true),
                                {})
                                )
                            )),
                        LForce(std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
                            *form_magForce_virtualWork_L, {P1}, {{"A", A}, {"B", B}, {"H", H}, {"nu", nuDiff}, {"Wmag0", Wmag0}, {"vacInd", vacInd}}, 
                                {   {"c0",std::make_shared<const fem::Constant<T>>(0.0)}, 
                                    {"c1",std::make_shared<const fem::Constant<T>>(0.0)}, 
                                    {"c2",std::make_shared<const fem::Constant<T>>(0.0)}}, {}, {}))),
                        domDepth(domDepthIn)
        {
            vacInd->x()->set(1.0);
            vacInd->name = "vacuumIndicator";
        }

    template<typename T> void force_calculation<T>::add_force_calc(const mag_tools::scen::xml_model_entity& fCalcParam){
        for (auto&i: determine_index_vector(fCalcParam)){
            this->forceCalculations.push_back(nodal_force<T>(i,this->F_vec, this->dofCoords, this->P1, this->boundaryMarkers, determine_rot_center(fCalcParam)));
        }
    }
    template<typename T> void force_calculation<T>::update_forces(const double& angle){
        assemble_force_vector();
        for (auto& fCalc: forceCalculations){
            fCalc.update_force(domDepth, angle);
        }
    }

    template<typename T> void force_calculation<T>::set_domain_depth(const double& domDepthIn){
        this->domDepth = domDepthIn;
    }

    template<typename T> void force_calculation<T>::store_results(const double& time){
        for (auto& fCalc: forceCalculations){
            fCalc.store_results(time);
        }
    }

    template<typename T> void force_calculation<T>::append_forces_to_output(mag_tools::output::result_xml& resXML){
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0){
            for (auto& fCalc: forceCalculations){
                auto fxNode = resXML.append_time_series("F", {"domain", "unit", "outputFormat"}, {std::to_string(fCalc.get_domain_idx()), "N", "xy"});
                fxNode->add_values(fCalc.get_force_results());
                auto TzNode = resXML.append_time_series("Tz", {"domain", "unit", "outputFormat"}, {std::to_string(fCalc.get_domain_idx()), "Nm", "z"});
                TzNode->add_values(fCalc.get_torque_results());
            }
        }
    }




    template<typename T> void force_calculation<T>::assemble_force_vector(){
        F_vec->set(0.0);
        dolfinx::fem::assemble_vector(F_vec->mutable_array(), *LForce);
    }
    

  

template class energy_density<PetscScalar>;
template class force_calculation<PetscScalar>;
template class nodal_force<PetscScalar>;
}

