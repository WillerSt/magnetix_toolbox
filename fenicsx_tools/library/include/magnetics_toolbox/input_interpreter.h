#pragma once

#include <iostream>

#include <dolfinx.h>
#include <dolfinx/io/XDMFFile.h>
#include <basix/e-lagrange.h>

#include <magnetics_toolbox/mag_tools_basic.h>
#include <magnetics_toolbox/coupled_material_model.h>
#include <magnetics_toolbox/field_sources.h>
#include <magnetics_toolbox/mesh_container.h>
#include <magnetics_toolbox/scenario_description.h>
#include <magnetics_toolbox/file_management.h>



#include <Eigen/Core>


namespace mag_tools{

    class solver_parameters{
        public:
        const double relTol;
        const double minRes;
        const int maxIter;

        public:
        solver_parameters(const mag_tools::scen::xml_model_entity& descIn);
    };

    class time_stepping{
        public:
        const std::vector<double> timeVec;

        private:
        const double* time;
        std::size_t curStep = 0;
        const std::size_t nStep;
        bool finished = false;

        public:

        time_stepping(const mag_tools::scen::xml_model_entity& descIn);

        void next_time_step();

        inline double get_time(){
            return *time;
        }
        inline std::size_t get_step(){
            return curStep;
        }
        inline std::size_t number_of_steps(){
            return this->nStep;
        }
        inline bool is_finished(){
            return finished;
        }

        inline std::vector<double> get_time_vec(){
            return this->timeVec;
        }

        double get_dt(){
            if (curStep > 0){
                return timeVec[curStep]-timeVec[curStep-1];
            }
            else{
                // initialization needs to be handeled outside of time stepping
                return 1.0;
            }
        }

        private:
        std::vector<double> construct_time_vec(const mag_tools::scen::xml_model_entity& descIn) const;  

    };

    class boundary_definition{
        public:
        const int idx;
        const double value;
        const std::string type;

        public: 
        boundary_definition(const mag_tools::scen::xml_model_entity& descIn);
    };

    
    template <typename T> class material_container{

        private:
        
        std::vector<std::shared_ptr<mag_tools::magnetic_material_coupling<T>>> materialModels = {};

        const int domainIndex;

        const std::shared_ptr<const mag_tools::quadrature_dof_coupler<T>> dofCouplerIn;
        const std::shared_ptr<const mag_tools::quadrature_dof_coupler<T>> dofCouplerOut;
        const std::shared_ptr<const mag_tools::quadrature_dof_coupler<T>> dofCouplerJac;

        std::string matType;

        bool energyCouplingInitialized = false;
        bool magneticVacuum = false;
        /*
        const auto dofCouplerM = std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(M, meshMarkers, idxIron, 2, 1);
        const auto dofCouplerNuDiff = std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(nuDiff, meshMarkers, idxIron, 2, 2);
        const auto dofCouplerB = std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(Bquad, meshMarkers, idxIron, 2, 1);
        */
       
        public:

        material_container(const mag_tools::scen::xml_model_entity& descIn, 
            const std::shared_ptr<dolfinxFunction<T>>& inputFunction, 
            const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
            const std::shared_ptr<dolfinxFunction<T>>& outputJac,
            const mesh_container<T>& meshContainer,
            const int& modelDim);

        material_container(const mag_tools::scen::xml_model_entity& descIn, 
            const std::shared_ptr<dolfinxFunction<T>>& inputFunction, 
            const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
            const std::shared_ptr<dolfinxFunction<T>>& outputJac,
            const std::shared_ptr<mesh_container<T>>& meshContainer,
            const int& modelDim);
        
        void update_material_models();

        void prepare_next_timestep();
        void set_diagonal(const T& diagVal);

        void couple_energy_density(const std::shared_ptr<const quadrature_dof_coupler<T>>& energyCoupler){
            for (size_t idx = 0; idx < materialModels.size(); idx++){
                materialModels[idx]->couple_energy_density((*energyCoupler).get_ref_vec());
            }
            energyCouplingInitialized = true;
        }

        void update_energy_density(){
            if (energyCouplingInitialized == true){
                for (auto& matModel: materialModels){
                    matModel->update_energy_density();
                }
            }
            else{
                std::cout << "WARNING: Requested energy density update but coupling was never initialized, no operation performed\n";
            }
        }

        std::string get_material_type(){
            return this->matType;
        }

        int get_domain_index() const{
            return domainIndex;
        }

        bool equiv_vacuum(){
            return magneticVacuum;
        }
        private:
        int determine_index(const mag_tools::scen::xml_model_entity& descIn);           

    };

    template <typename T> class source_container{
    

        std::shared_ptr<mag_tools::src::const_scalar_source<T>> source;



        public:
        source_container(const mag_tools::scen::xml_model_entity& descIn, 
            const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
            const mesh_container<T>& meshContainer);

        void update_source(const double& t);

        inline std::shared_ptr<mag_tools::src::const_scalar_source<T>> get_source(){
            return this->source;
        }

        

    };

    template <typename T> std::vector<int> get_magnetic_domains(std::shared_ptr<std::vector<mag_tools::material_container<T>>> matDefs){
        std::vector<int> vacIdx;
        for (auto& mat: *matDefs){
            if(!mat.equiv_vacuum()){
                vacIdx.push_back(mat.get_domain_index());
            }            
        }
        std::sort(vacIdx.begin(), vacIdx.end());
        vacIdx.erase(std::unique(vacIdx.begin(), vacIdx.end()), vacIdx.end());
        return vacIdx;
    }
}
