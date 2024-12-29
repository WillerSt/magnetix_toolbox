#pragma once
#include <magnetics_toolbox/input_interpreter.h>
#include <magnetics_toolbox/output.h>

namespace mag_tools{


template<typename T> class nodal_force{
    //using T = PetscScalar;
    using EigenVec = Eigen::Matrix<double, 2, 1>;

    const std::size_t idx;
    const std::vector<int32_t> fDofs;
    const std::shared_ptr<dolfinxRHS<T>> fVec;

    std::vector<EigenVec> nodalForces;
    std::vector<EigenVec> nodeCoords;

    EigenVec rotCenter;
    EigenVec F = {0.0, 0.0};
    double Tz = 0.0;

    std::vector<EigenVec> forceSeries;
    std::vector<double> torqueSeries;
    std::vector<double> timeSeries;

    private:
        const std::shared_ptr<std::vector<double>> dofCoords;

    public:
        nodal_force(const int& bdryIdx, const std::shared_ptr<dolfinxRHS<T>>& fVecIn,
                const std::shared_ptr<std::vector<double>>& dofCoords,
                const std::shared_ptr<const dolfinxFS<T>>& P1, 
                const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>>& boundaryMarkers,
                const EigenVec& rotCenterIn = {0,0});

        void update_force(const double& scaleFac, const double& angle);

        void store_results(const double& time);

        size_t get_domain_idx() const{
            return this->idx;
        }

        std::vector<double> get_torque_results() const{
            return torqueSeries;
        }
        std::vector<EigenVec> get_force_results() const{
            return forceSeries;
        }
};


template<typename T> class energy_density{
    // using T = PetscScalar;

    private:
    const std::shared_ptr<std::vector<mag_tools::material_container<T>>> matDefs;

    const std::shared_ptr<dolfinxMesh<T>> mesh;
    const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> meshMarkers;

    const std::shared_ptr<dolfinxFunction<T>> Wmag;

    std::vector<std::shared_ptr<mag_tools::quadrature_dof_coupler<T>>> dofCouplers;

    const std::shared_ptr<dolfinxFunction<T>> B;
    const std::shared_ptr<dolfinxFunction<T>> H;

    const quadrature_dof_coupler<T> WmagCplg = quadrature_dof_coupler<T>(Wmag,1,1);
    const quadrature_dof_coupler<T> BCplg = quadrature_dof_coupler<T>(B,2,1);
    const quadrature_dof_coupler<T> HCplg = quadrature_dof_coupler<T>(H,2,1);   

    const std::vector<std::vector<std::vector<double*>>> WmagDofs = *WmagCplg.get_ref_vec() ; 
    const std::vector<std::vector<std::vector<double*>>> BDofs = *HCplg.get_ref_vec() ; 
    const std::vector<std::vector<std::vector<double*>>> HDofs = *BCplg.get_ref_vec() ; 

    public:

    energy_density(const mesh_container<T>& meshContainer,
                    const std::shared_ptr<dolfinxFunction<T>>& BIn,
                    const std::shared_ptr<dolfinxFunction<T>>& HIn,
                    const std::shared_ptr<std::vector<mag_tools::material_container<T>>>matDefsIn );
    void perform_update();

    std::shared_ptr<dolfinxFunction<T>> get_density_function(){
        return Wmag;
    }    
};

template<typename T> class force_calculation{
    private:
    //using T = PetscScalar;
    const std::shared_ptr<dolfinxFunction<T>> A;
    const std::shared_ptr<dolfinxFunction<T>> B;
    const std::shared_ptr<dolfinxFunction<T>> H;
    const std::shared_ptr<dolfinxFunction<T>> nuDiff;

    const std::shared_ptr<dolfinxMesh<T>> mesh;
     const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> meshMarkers;
    const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> boundaryMarkers;

    const std::shared_ptr<const dolfinxFS<T>> P1;

    const std::shared_ptr<dolfinxFunction<T>> Wmag0;

    const std::shared_ptr<dolfinxFunction<T>> vacInd;

    const std::shared_ptr<const fem::Form<T>> LForce;

    const std::shared_ptr<dolfinxRHS<T>> F_vec = std::make_shared<dolfinxRHS<T>>(dolfinx::la::Vector<T>(LForce->function_spaces()[0]->dofmap()->index_map,
                LForce->function_spaces()[0]->dofmap()->index_map_bs()));

    const std::shared_ptr<std::vector<double>>  dofCoords = std::make_shared<std::vector<double>>(P1->tabulate_dof_coordinates(false));

    std::vector<nodal_force<T>> forceCalculations;

   double domDepth;


    public:
    force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                        const std::shared_ptr<dolfinxFunction<T>>& BIn,
                        const std::shared_ptr<dolfinxFunction<T>>& HIn,
                        const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                        const std::shared_ptr<energy_density<T>>& energyDensity,
                        const mag_tools::mesh_container<T>& meshContainer,
                        const std::vector<std::size_t>& idxIn,
                        const double& domDepthIn,
                        Eigen::Matrix<double, 2, 1> rotCenter= {0.0,  0.0});

        force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                        const std::shared_ptr<dolfinxFunction<T>>& BIn,
                        const std::shared_ptr<dolfinxFunction<T>>& HIn,
                        const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                        const std::shared_ptr<energy_density<T>>& energyDensity,
                        const mag_tools::mesh_container<T>& meshContainer,
                        const mag_tools::scen::xml_model_entity& fCalcParam);

        force_calculation(const std::shared_ptr<dolfinxFunction<T>>& AIn,
                    const std::shared_ptr<dolfinxFunction<T>>& BIn,
                    const std::shared_ptr<dolfinxFunction<T>>& HIn,
                    const std::shared_ptr<dolfinxFunction<T>>& nuDiffIn,
                    const std::shared_ptr<energy_density<T>>& energyDensity,
                    const mag_tools::mesh_container<T>& meshContainer,
                    const double& domDepthIn = 1.0);

    void add_force_calc(const mag_tools::scen::xml_model_entity& fCalcParam);
    void update_forces(const double& angle = 0.0);
    void set_domain_depth(const double& domDepthIn);
    void set_mag_domains(const std::vector<int>& magDoms){
        std::cout << "mag domains:";
        for (auto& idx: magDoms){
            // defined for quadrature function space, stencil is dg0
            auto vacSetter = quadrature_dof_coupler<T>(vacInd, this->meshMarkers, idx, 1, 1);
            vacSetter.set_all_entries(0.0);
            std::cout << " " << idx;
        }
        std::cout << std::endl;
    }
    void store_results(const double& time);

    void append_forces_to_output(mag_tools::output::result_xml& resXML);
    std::shared_ptr<dolfinxFunction<T>> get_w0_function(){
        return Wmag0;
    }

    void output_indicator_function(const std::string& fileName){
        io::VTKFile file(MPI_COMM_WORLD, fileName+".pvd", "w");
        file.write<T>({*vacInd}, 0.0);
        
        //auto outFile = std::make_unique<dolfinx::io::VTXWriter<T>>(MPI_COMM_WORLD, fileName, dolfinx::io::adios2_writer::U<T>({this->vacInd}));
        //outFile->write(0.0);
        io::VTXWriter<U<T>> outFile(MPI_COMM_WORLD, fileName+".bp", {vacInd}, "bp4");
        outFile.write(0.0);
        
    }
    private:
    void assemble_force_vector();
    
    public:
    static std::vector<std::size_t> determine_index_vector(const mag_tools::scen::xml_model_entity& params){
        auto domainVec = scen::search_for_vector(params.get_parameters(), "domain");
        std::vector<std::size_t> idxVec;
        for (auto& dom: domainVec){
            idxVec.push_back(scen::as_int(params, dom));
        }
        return idxVec;
    }
    static Eigen::Matrix<double, 2, 1> determine_rot_center(const mag_tools::scen::xml_model_entity& params){
        Eigen::Matrix<double, 2, 1> rotCenter = {scen::as_double(params, "rotCenter_x"), scen::as_double(params, "rotCenter_y")};
        return rotCenter;
    }



        

};



}

