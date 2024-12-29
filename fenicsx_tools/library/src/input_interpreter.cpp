#include <magnetics_toolbox/input_interpreter.h>
#include <magnetics_toolbox/mesh_container.h>
#include <magnetics_toolbox/domain_movement.h>

namespace mag_tools{

    template<typename T> mesh_container<T>::mesh_container(const mag_tools::scenario_description& desc):
        /*
        coordElement(   dolfinx::fem::CoordinateElement<U<T>>(
                            std::make_shared<basix::FiniteElement<T>>(
                                basix::element::create_lagrange<T>(
                                    basix::cell::type::triangle, 1,basix::element::lagrange_variant::unset, false)))),*/
        coordElement( dolfinx::fem::CoordinateElement<U<T>>(dolfinx::mesh::CellType::triangle,1, basix::element::lagrange_variant::unset)),
        /*coordElement(std::make_shared<basix::FiniteElement<T>>(
          basix::create_element<T>(basix::element::family::P,
                                   dolfinx::mesh::cell_type_to_basix_type(dolfinx::mesh::CellType::triangle),
                                   1, basix::element::lagrange_variant::unset,
                                  basix::element::dpc_variant::unset, false))),*/
        meshInput(dolfinx::io::XDMFFile(MPI_COMM_WORLD, desc.meshInput.mesh_file(), "r")),
        boundaryInput(dolfinx::io::XDMFFile(MPI_COMM_WORLD, desc.meshInput.boundary_file(), "r")),
        mesh(this->initialize_mesh()),
        meshMarkers(std::make_shared<const dolfinx::mesh::MeshTags<std::int32_t>>(this->meshInput.read_meshtags(*mesh,"Grid"))),
        boundaryMarkers(std::make_shared<const dolfinx::mesh::MeshTags<std::int32_t>>(this->boundaryInput.read_meshtags(*mesh,"Grid"))){

    }

    template<typename T> std::shared_ptr<dolfinxMesh<T>> mesh_container<T>::initialize_mesh(){
        std::shared_ptr<dolfinxMesh<T>> meshTemp = std::make_shared<dolfinxMesh<T>>(this->meshInput.read_mesh(this->coordElement, dolfinx::mesh::GhostMode::shared_facet, "Grid"));
        // needed to identify boundaries
        meshTemp->topology()->create_connectivity(mesh->topology()->dim()-1,mesh->topology()->dim());
        // needed for coordinate shifts
        meshTemp->topology()->create_connectivity(mesh->topology()->dim(), 0);
        return meshTemp;
    }



    solver_parameters::solver_parameters(const mag_tools::scen::xml_model_entity& descIn):
        relTol(std::stod(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(),"rTol")].value.c_str())), 
        minRes(std::stod(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(),"minRes")].value.c_str())), 
        maxIter(std::stoi(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(),"maxIter")].value.c_str())){

    }


    time_stepping::time_stepping(const mag_tools::scen::xml_model_entity& descIn):
        timeVec(construct_time_vec(descIn)),
        time(&timeVec[0]),
        nStep(timeVec.size()){

    }

    std::vector<double> time_stepping::construct_time_vec(const mag_tools::scen::xml_model_entity& descIn) const{
        std::vector<double> timeSteps;

        auto parameters =  descIn.get_parameters();
        std::string type = parameters[scen::search_for_name(parameters, "type")].value;

        if (strcmp(type.c_str(), "linear")==0){
            double timeStart = std::stod(parameters[scen::search_for_name(parameters, "tStart")].value);
            double timeEnd = std::stod(parameters[scen::search_for_name(parameters, "tEnd")].value);
            int nSteps = std::stoi(parameters[scen::search_for_name(parameters, "nSteps")].value);
            double timeStep = (timeEnd - timeStart)/(nSteps-1);
            
            if (!(timeStep > 0))
            {
                std::cout << "Time step <= 0: " << timeStep << std::endl; 
            }

            for (int i = 0; i < nSteps; i++){
                timeSteps.push_back(timeStart + i*timeStep);
            }

        }
        else if (strcmp(type.c_str(), "xmlFile")==0){
            auto seriesXML = parameters[scen::search_for_name(parameters, "xmlFile")].value;
            
            if (!boost::filesystem::exists(seriesXML)){
                std::cout << "ERROR: Time series XML " << seriesXML << " does not exist, aborting\n";
                exit(1);
            }
            else{
                timeSteps = scen::read_xml_series(seriesXML, "time");
                /*
                std::cout << "Read timesteps:";
                for (std::size_t i = 0; i < timeSteps.size(); i++){
                    std::cout <<" " << timeSteps[i];
                }
                std::cout << std::endl;
                */
            }
        }
        else if (strcmp(type.c_str(), "singleStep")==0){
            timeSteps = {std::stod(parameters[scen::search_for_name(parameters, "tStart")].value)};
        }
        else{
            std::cout << "Time stepping " << type << " not implemented\n";
            exit(1); 
        }




        return timeSteps;
    }

    void time_stepping::next_time_step(){
        if (curStep < (nStep-1)){
            curStep ++;
            time = &timeVec[curStep];
        }
        else{
            finished = true;
        }
    }

    boundary_definition::boundary_definition(const mag_tools::scen::xml_model_entity& descIn):
        idx(std::stoi(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(), "index")].value.c_str())),
        value(std::stod(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(), "value")].value.c_str())),
        type(descIn.get_parameters()[scen::search_for_name(descIn.get_parameters(), "type")].value){

    }

    template <typename T> material_container<T>::material_container(const mag_tools::scen::xml_model_entity& descIn, 
        const std::shared_ptr<dolfinxFunction<T>>& inputFunction, 
        const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
        const std::shared_ptr<dolfinxFunction<T>>& outputJac,
        const mesh_container<T>& meshContainer,
        const int& modelDim):domainIndex(determine_index(descIn)),
        dofCouplerIn(std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(inputFunction, meshContainer.get_mesh_markers(), domainIndex, modelDim, 1)),
        dofCouplerOut(std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(outputFunction, meshContainer.get_mesh_markers(), domainIndex, modelDim, 1)),
        dofCouplerJac(std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(outputJac, meshContainer.get_mesh_markers(), domainIndex, modelDim, modelDim)){
        
        auto parameters =  descIn.get_parameters();
        std::string type = parameters[scen::search_for_name(parameters, "type")].value;
        this->matType = type;

        // hysteresis model    
        if ((strcmp(type.c_str(), "HGM") == 0)|| (strcmp(type.c_str(), "hgm")==0) ){

            std::string dpcXML = parameters[scen::search_for_name(parameters, "distFile")].value;
            file_manager::abort_if_not_existent(dpcXML, "material_container");

            double Bsat = std::stod(parameters[scen::search_for_name(parameters, "Bsat")].value);

            auto gridConstructor = xml_grid_constructor<2, CS_sphere_cont<2>>(dpcXML);
            const auto mallocGrid = std::make_shared<const malloc_grid<2, CS_sphere_cont<2>>>(gridConstructor, Bsat);
            auto dpcParameter = mag_tools::parameters_DPC_Model_Cont_2D(mallocGrid);
            
            for (std::size_t i = 0; i<dofCouplerOut->get_ref_vec()->size(); i++){
                materialModels.push_back(std::make_shared<mag_tools::coupled_DPC_Model_Cont_2D>(dpcParameter, dofCouplerIn->get_ref_vec(), dofCouplerJac->get_ref_vec(), dofCouplerOut->get_ref_vec(), i));
            }

        } // atan model
        else if((strcmp(type.c_str(), "ATAN")==0) || (strcmp(type.c_str(), "atan")==0)){
            auto atanParams =  mag_tools::parameters_atan_H_B(
                std::stod(parameters[scen::search_for_name(parameters, "Bsat")].value), 
                std::stod(parameters[scen::search_for_name(parameters, "muInit")].value),
                std::stod(parameters[scen::search_for_name(parameters, "muLin")].value));

            for (std::size_t i = 0; i<dofCouplerOut->get_ref_vec()->size(); i++){
                materialModels.push_back(std::make_shared<mag_tools::coupled_atan_H_B_2D>(atanParams, dofCouplerIn->get_ref_vec(), dofCouplerJac->get_ref_vec(), dofCouplerOut->get_ref_vec(), i));
            }
            
        } // linear material
        else if (strcmp(type.c_str(), "linear")==0){
            auto linParams = mag_tools::parameters_lin_H_B(std::stod(parameters[scen::search_for_name(parameters, "muR")].value));

            for (std::size_t i = 0; i<dofCouplerOut->get_ref_vec()->size(); i++){
                materialModels.push_back(std::make_shared<mag_tools::coupled_lin_H_B>(linParams, dofCouplerIn->get_ref_vec(), dofCouplerJac->get_ref_vec(), dofCouplerOut->get_ref_vec(), i));
            }
        } // permanent magnet
        else if (strcmp(type.c_str(), "permanentMagnet")==0){
            auto model = parameters[scen::search_for_name(parameters, "model")].value.c_str();
            if (strcmp(model, "linear")==0){

                auto Br = std::stod(parameters[scen::search_for_name(parameters, "Br")].value.c_str());
                auto muR = std::stod(parameters[scen::search_for_name(parameters, "muR")].value.c_str());
                auto xDir = std::stod(parameters[scen::search_for_name(parameters, "xDir")].value.c_str());
                auto yDir = std::stod(parameters[scen::search_for_name(parameters, "yDir")].value.c_str());
                
                auto magParams = mag_tools::parameters_lin_PM_2D(muR, Br, xDir, yDir);
                for (std::size_t i = 0; i<dofCouplerOut->get_ref_vec()->size(); i++){
                    materialModels.push_back(std::make_shared<mag_tools::coupled_const_mag>(magParams, dofCouplerIn->get_ref_vec(), dofCouplerJac->get_ref_vec(), dofCouplerOut->get_ref_vec(), i));
                }
            }            
        }
        else if (strcmp(type.c_str(), "table")==0){
            auto xmlFile = mag_tools::scen::as_string(descIn, "xmlFile");
            auto Htable = scen::read_xml_material(xmlFile, "H");
            auto Btable = scen::read_xml_material(xmlFile, "B");
            auto splineParameters = mag_tools::parameters_spline_H_B(Htable, Btable);
            for (std::size_t i = 0; i<dofCouplerOut->get_ref_vec()->size(); i++){
                materialModels.push_back(std::make_shared<mag_tools::coupled_spline_H_B_2D>(splineParameters, dofCouplerIn->get_ref_vec(), dofCouplerJac->get_ref_vec(), dofCouplerOut->get_ref_vec(), i));
            }
        }

    }

    template<typename T> material_container<T>:: material_container(const mag_tools::scen::xml_model_entity& descIn, 
        const std::shared_ptr<dolfinxFunction<T>>& inputFunction, 
        const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
        const std::shared_ptr<dolfinxFunction<T>>& outputJac,
        const std::shared_ptr<mesh_container<T>>& meshContainer,
        const int& modelDim):material_container(descIn, inputFunction, outputFunction, outputJac, *meshContainer, modelDim){

    }
    
    template<typename T> void material_container<T>::update_material_models(){
        for (auto& matCplng: materialModels){
            matCplng->update_coupling_entries_M_nu();
        }
    }

    template<typename T> void material_container<T>::prepare_next_timestep(){
        for (auto& matCplng: materialModels){
            matCplng->accept_state();
        }
    }
    
    template<typename T> void material_container<T>::set_diagonal(const T& diagVal){
        dofCouplerJac->set_diagonal(diagVal);
    }

    template<typename T> int material_container<T>::determine_index(const mag_tools::scen::xml_model_entity& descIn){
        auto parameters =  descIn.get_parameters();
        for (std::size_t i = 0; i<parameters.size(); i ++){
            if (std::strcmp("domainIndex", parameters[i].name.c_str()) == 0){
                int domainIndex = stoi(parameters[i].value);
                return domainIndex;
            }
            
        }
        std::cout << "Error: No domain index given \n";
        exit(1);
    }  

    template<typename T> source_container<T>::source_container(
            const mag_tools::scen::xml_model_entity& descIn, 
            const std::shared_ptr<dolfinxFunction<T>>& outputFunction,
            const mesh_container<T>& meshContainer){
            

        auto parameters =  descIn.get_parameters();
        std::string type = parameters[scen::search_for_name(parameters, "type")].value;

        if (strcmp(type.c_str(), "currentDensity") ==0){
            if (strcmp(parameters[scen::search_for_name(parameters, "distribution")].value.c_str(), "const") == 0){
                auto evolution = parameters[scen::search_for_name(parameters, "evolution")].value.c_str();
                if (strcmp(evolution, "sine") == 0){
                    auto domainVec  = search_for_vector(parameters, "domain");
                    auto scaleVec  = search_for_vector(parameters, "scale");

                    std::vector<double> scale;
                    if ((domainVec.size() < 1) || (scaleVec.size()<1) || (scaleVec.size()!=domainVec.size())){
                        std::cout << "Error in source definition scale and domain information inconsistent" << std::endl;
                        exit(1);
                    }
                    for (std::size_t i = 0; i < scaleVec.size(); i++){
                        scale.push_back(std::stod(parameters[scen::search_for_name(parameters, scaleVec[i])].value.c_str()));
                    }
                    std::vector<int> domain;
                    for (std::size_t i = 0; i < domainVec.size(); i++){
                        domain.push_back(std::stoi(parameters[scen::search_for_name(parameters, domainVec[i])].value.c_str()));
                    }


                    auto sineAmp = std::stod(parameters[scen::search_for_name(parameters, "sineAmp")].value.c_str());
                    auto sinePhase = std::stod(parameters[scen::search_for_name(parameters, "sinePhase")].value.c_str());
                    auto sineOffset = std::stod(parameters[scen::search_for_name(parameters, "sineOffset")].value.c_str());
                    auto sineFreq = std::stod(parameters[scen::search_for_name(parameters, "sineFreq")].value.c_str());

                    source = std::make_shared<mag_tools::src::const_scalar_source<T>>(outputFunction, meshContainer.get_mesh_markers(), scale, domain,
                            std::make_shared<const mag_tools::src::sine_variation>(sineAmp, sineFreq, sinePhase, sineOffset));
                    
                }
                else if(strcmp(evolution, "timeSeries")==0){
                    auto domainVec  = search_for_vector(parameters, "domain");
                    auto scaleVec  = search_for_vector(parameters, "scale");

                    std::vector<double> scale;
                    if ((domainVec.size() < 1) || (scaleVec.size()<1) || (scaleVec.size()!=domainVec.size())){
                        std::cout << "Error in source definition scale and domain information inconsistent" << std::endl;
                        exit(1);
                    }
                    for (std::size_t i = 0; i < scaleVec.size(); i++){
                        scale.push_back(std::stod(parameters[scen::search_for_name(parameters, scaleVec[i])].value.c_str()));
                    }
                    std::vector<int> domain;
                    for (std::size_t i = 0; i < domainVec.size(); i++){
                        domain.push_back(std::stoi(parameters[scen::search_for_name(parameters, domainVec[i])].value.c_str()));
                    }

                    auto seriesXML = parameters[scen::search_for_name(parameters, "xmlFile")].value;
                    auto quantity = parameters[scen::search_for_name(parameters, "quantityName")].value;

                    auto timeSteps = scen::read_xml_series(seriesXML, "time");
                    /*
                    std::cout << "Read timesteps:";
                    for (std::size_t i = 0; i < timeSteps.size(); i++){
                        std::cout <<" " << timeSteps[i];
                    }
                    */
                    auto inputValues = scen::read_xml_series(seriesXML, quantity.c_str());             
                    /*
                    std::cout << "\nRead inut:";
                    for (std::size_t i = 0; i < inputValues.size(); i++){
                        std::cout <<" " << inputValues[i];
                    }
                    std::cout << std::endl;
                    */
                    source = std::make_shared<mag_tools::src::const_scalar_source<T>>(outputFunction, meshContainer.get_mesh_markers(), scale, domain,
                            std::make_shared<const mag_tools::src::value_table>(timeSteps, inputValues));
                    

                     


                }
                else{
                    std::cout << "ERROR: Source evolution " << evolution << " is not implemented" << std::endl;
                    exit(1);
                }
            }
        }

    }

    template<typename T> void source_container<T>::update_source(const double& t){
        this->source->update_source(t);
    }

    template class mesh_container<PetscScalar>;
    template class material_container<PetscScalar>;
    template class source_container<PetscScalar>;
}