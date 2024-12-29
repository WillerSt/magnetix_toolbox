#include <magnetics_toolbox/maxwell_solvers.h>
#include <magnetics_toolbox/mag_tools_basic.h>
#include <dolfinx/fem/Constant.h>

#include<dolfinx/io/XDMFFile.h>
#include<dolfinx/io/VTKFile.h>
#include<dolfinx/io/ADIOS2Writers.h>

#include <A_form_NR.h>


#include <iostream>
#include <math.h>

#include <magnetics_toolbox/linear_solver.h>

#include <magnetics_toolbox/coupled_material_model.h>

#include <magnetics_toolbox/field_sources.h>
#include <magnetics_toolbox/field_quantity.h>

#include <magnetics_toolbox/input_interpreter.h>

#include <magnetics_toolbox/function_processing.h>

#include <magnetics_toolbox/file_management.h>
#include <boost/program_options.hpp>

#include <magnetics_toolbox/post_processing.h>
#include <magForce_virtualWork.h>

#include <magnetics_toolbox/utilities.h>
#include <magnetics_toolbox/utility.h>
#include <basix/quadrature.h>
#include <basix/e-lagrange.h>

namespace mag_tools{
    void magnetic_field_2D_Az(int argc, char* argv[]){
        using T = PetscScalar;

        dolfinx::init_logging(argc,  argv);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        auto scenFile = mag_tools::get_scenario_file(argc, argv);


        auto scen = mag_tools::scenario_description(scenFile);

        mag_tools::file_manager fileManager(scen.get_scen_path(), scen.scenario_name(), true);
        
        // model parameter
        const double nu0 = constants::mu0inv;

        // config
        
        auto meshContainer = mag_tools::mesh_container<T>(scen);
        fileManager.copy_all_input_files(scen);
        
        const auto mesh = meshContainer.get_mesh();
        const auto meshMarkers = meshContainer.get_mesh_markers();
        const auto boundaryMarkers = meshContainer.get_boundary_markers();

        const std::string fieldsFile = (fileManager.result_path()/"fields_BH.bp").c_str();
        const std::string potentialFile = (fileManager.result_path()/"potential_Az.bp").c_str();
        const std::string currDensFile = (fileManager.result_path()/"currDens_j.bp").c_str();
        const std::string magFile = (fileManager.result_path()/"fields_M.bp").c_str();

        // initialize function spaces and functions
        const size_t dimG = mesh->topology()->dim();

        /*
        auto quadRule = basix::quadrature::make_quadrature<T>(basix::quadrature::type::Default, basix::cell::type::triangle, basix::polyset::type::standard, degQ);

        auto qScaE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, 1));
        
        auto qVecE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, 2));
        auto qMatE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, 4));
        
        //auto qVecE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>({qScaE ,qScaE }));
        //auto qMatE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>({qScaE, qScaE, qScaE, qScaE }));
        
        auto const dofmapSca = std::make_shared<const dolfinx::fem::DofMap>(
                                    dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                                    dolfinx::fem::create_element_dof_layout<T>(*qScaE), 
                                    *mesh->topology(), nullptr, nullptr));

        auto const dofmapVec = std::make_shared<const dolfinx::fem::DofMap>(
                            dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                            dolfinx::fem::create_element_dof_layout<T>(*qVecE), 
                            *mesh->topology(), nullptr, nullptr));

        auto const dofmapMat = std::make_shared<const dolfinx::fem::DofMap>(
                            dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                            dolfinx::fem::create_element_dof_layout<T>(*qMatE), 
                            *mesh->topology(), nullptr, nullptr));

        const auto quadScaFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
            dolfinx::fem::FunctionSpace<T>(mesh, qScaE, dofmapSca,{1}));

        const auto quadVecFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
            dolfinx::fem::FunctionSpace<T>(mesh, qVecE, dofmapVec,{2}));

        const auto quadMatFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
            dolfinx::fem::FunctionSpace<T>(mesh, qMatE, dofmapMat,{2,2}));

    */
        const auto quadScaFS = utility::create_quadrature_functionspace<T>(mesh, "sca");
        const auto quadVecFS = utility::create_quadrature_functionspace<T>(mesh, "vec");
        const auto quadMatFS = utility::create_quadrature_functionspace<T>(mesh, "mat");

        
        auto quadScaFunc = std::make_shared<fem::Function<T>>(quadScaFS);
        auto quadVecFunc = std::make_shared<fem::Function<T>>(quadVecFS);
        auto quadMatFunc = std::make_shared<fem::Function<T>>(quadMatFS);
        
        auto dg1Elemem = basix::create_element<U<T>>(
            basix::element::family::P, basix::cell::type::triangle, 1,
            basix::element::lagrange_variant::unset,
            basix::element::dpc_variant::unset, true);        

        const auto dgVecFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {dimG,1}));
        const auto dgScaFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {}));      

        auto cg2Elemem = basix::create_element<U<T>>(
        basix::element::family::P, basix::cell::type::triangle, 2,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, false);

        auto cgScaFS = std::make_shared<fem::FunctionSpace<T>>(fem::create_functionspace(mesh, cg2Elemem, {}));
        
        const auto nuDiff = std::make_shared<dolfinxFunction<T>>(quadMatFS);     
        const auto M  = std::make_shared<dolfinxFunction<T>>(quadVecFS); 
        const auto M0 = std::make_shared<dolfinxFunction<T>>(quadVecFS);     
        const auto j  = std::make_shared<dolfinxFunction<T>>(quadScaFS); 
        const auto A  = std::make_shared<dolfinxFunction<T>>(cgScaFS);
        const auto A_delta = std::make_shared<dolfinxFunction<T>>(cgScaFS);

        const auto zeroCG = mag_tools::create_const_function<T>(cgScaFS,T(0.0));
        const auto zeroQuad  = mag_tools::create_const_function<T>(quadVecFS, T(0.0));

        M->x()->set(0.0);
        nuDiff->x()->set(0.0);
        j->x()->set(0.0);
        A_delta->x()->set(0.0);
        A->x()->set(0.0);
        M0->x()->set(0.0);

        auto _A_delta = dolfinxVector(la::petsc::create_vector_wrap(*A_delta->x()), false);
        auto _A = dolfinxVector(la::petsc::create_vector_wrap(*A->x()), false);
        
        // evaluation of B
        const std::shared_ptr<mag_tools::field_quantity<T>> curlEval = std::make_shared<mag_tools::surf_curl_evaluation<T>>(A, quadVecFS, 1.0);
        const auto Bquad = curlEval->get_result_function();

        // field quantities
        std::shared_ptr<mag_tools::field_quantity<T>> magFieldCalc = std::make_shared<mag_tools::linear_superposition<T>>(
            std::vector<std::shared_ptr<dolfinxFunction<T>>>({Bquad, M}), std::vector<double>({nu0, -1.0}));
        const auto Hquad = magFieldCalc->get_result_function();

        // output of B and H
        auto B = std::make_shared<dolfinxFunction<T>>(dgVecFS);
        auto H = std::make_shared<dolfinxFunction<T>>(dgVecFS);
        auto M_out = std::make_shared<dolfinxFunction<T>>(dgVecFS);
        auto M0_out = std::make_shared<dolfinxFunction<T>>(dgVecFS);


        auto j_out = std::make_shared<dolfinxFunction<T>>(dgScaFS);

        auto dgProj = mag_tools::dg_projection<T>(mesh, M->function_space()->mesh()->topology()->dim());
        auto dgProjSca = mag_tools::dg_projection<T>(mesh, 1);


        // renaming of function (has to be done before passing to output-file)
        nuDiff->name = "nuDiff";
        M->name = "M";
        M0->name = "M0";
        j->name = "j";
        A_delta->name = "A_delta";
        A->name = "A";
        Bquad->name = "B_quad";

        B->name = "B";
        H->name = "H";
        j_out->name = "j";
        M_out->name = "M";
        M0_out->name = "M0";

        // TODO:: Remove print functions __DEBUG
        //std::cout << "Sub elements vec:" << quadVecFS->element()->num_sub_elements() << " Mat element" << quadMatFS->element()->num_sub_elements() << std::endl;
        //std::cout<< "Cell: "<<  mesh->topology()->index_map(dimG)->size_global() <<" Dimension check vec: "<< M->x()->array().size() << " mat: " << nuDiff->x()->array().size() << std::endl;  
        // initialization of output file

        std::vector<std::unique_ptr<dolfinx::io::VTXWriter<T>>> outputFiles;
        outputFiles.push_back(std::make_unique<dolfinx::io::VTXWriter<T>>(MPI_COMM_WORLD, fieldsFile, dolfinx::io::adios2_writer::U<T>({B,H})));
        outputFiles.push_back(std::make_unique<dolfinx::io::VTXWriter<T>>(MPI_COMM_WORLD, magFile, dolfinx::io::adios2_writer::U<T>({M_out, M0_out})));
        outputFiles.push_back(std::make_unique<dolfinx::io::VTXWriter<T>>(MPI_COMM_WORLD, potentialFile, dolfinx::io::adios2_writer::U<T>({A})));
        outputFiles.push_back(std::make_unique<dolfinx::io::VTXWriter<T>>(MPI_COMM_WORLD, currDensFile, dolfinx::io::adios2_writer::U<T>({j_out})));
    
        // coupling of global jacobian
        const auto dofCouplerNuDiffAll = std::make_shared<const mag_tools::quadrature_dof_coupler<T>>(
                nuDiff, nuDiff->function_space()->mesh()->topology()->dim(), nuDiff->function_space()->mesh()->topology()->dim());

        // set permeability
        dofCouplerNuDiffAll->set_all_entries_zero();
        dofCouplerNuDiffAll->set_diagonal(nu0);

        // material setup, constant permanent magnet magnetization is stored in M0
        std::shared_ptr<std::vector<mag_tools::material_container<T>>> materialDefinition = 
            std::make_shared<std::vector<mag_tools::material_container<T>>>(); 
        for(auto& src: scen.materialInput){
            if (strcmp(src.get_type_parameter().c_str(), "permanentMagnet")==0){
                materialDefinition->push_back(mag_tools::material_container<T>(src, Bquad, M0, nuDiff, meshContainer, Bquad->function_space()->mesh()->topology()->dim()));
            }
            else{
                materialDefinition->push_back(mag_tools::material_container<T>(src, Bquad, M, nuDiff, meshContainer, Bquad->function_space()->mesh()->topology()->dim()));
            }
        }

        // source setup, current sources
        std::vector<mag_tools::source_container<T>> fieldSources;
        for(auto& src: scen.sourceInput){
            fieldSources.push_back(mag_tools::source_container<T>(src, j, meshContainer));
        }

        // initialize constants
        const auto nuConst = std::make_shared<const fem::Constant<T>>(nu0);

        // setup forms and corresponding linear system
        const auto a = std::make_shared<const varForm<T>>(dolfinx::fem::create_form<T>(
            *form_A_form_NR_a, {cgScaFS, cgScaFS}, {{"J", nuDiff}}, {}, {}, {}));
        const auto L = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
            *form_A_form_NR_L, {cgScaFS}, {{"j", j}, {"M0", M0}, {"M", M}, {"A_old", A}}, {{"nu", nuConst}}, {}, {}));

        const auto zeroL = std::make_shared<const fem::Form<T>>(dolfinx::fem::create_form<T>(
            *form_A_form_NR_L, {cgScaFS}, {{"j", j}, {"M0", M0}, {"M", zeroQuad}, {"A_old", zeroCG}}, {{"nu", nuConst}}, {}, {}));

        // setup boundaries
        std::vector<std::shared_ptr<const fem::DirichletBC<T>>> bc;
        for(auto& bdryDesc: scen.boundaryDef){
            mag_tools::boundary_definition bdryDef(bdryDesc);
            if (std::strcmp(bdryDef.type.c_str(), "Dirichlet") == 0){
                const auto bdofs = fem::locate_dofs_topological(*(cgScaFS->mesh()->topology_mutable()),
                        *(cgScaFS->dofmap()), cgScaFS->mesh()->topology()->dim()-1, boundaryMarkers->find(bdryDef.idx));
                bc.push_back(std::make_shared<const fem::DirichletBC<T>>(bdryDef.value, bdofs, cgScaFS));
            }
            
        }

        // setup energy calculation
        std::shared_ptr<energy_density<T>> wmag = std::make_shared<energy_density<T>>(meshContainer,Bquad, Hquad, materialDefinition);
        
        // setup force calculation
        std::vector<std::shared_ptr<force_calculation<T>>> forceCalculation;
        auto magDom = get_magnetic_domains<T>(materialDefinition);
        size_t cnt = 0;
        if (scen.forceCalc.size() > 0){
            auto Fcalc = std::make_shared<force_calculation<T>>(force_calculation<T>(A, Bquad, Hquad, nuDiff, wmag, meshContainer));
            forceCalculation.push_back(Fcalc);
            Fcalc->set_domain_depth(scen::as_double(scen.forceCalc[0], "domainDepth"));
            Fcalc->set_mag_domains(magDom);
            for(auto& FcalcDef: scen.forceCalc){
                Fcalc->add_force_calc(FcalcDef);
            }
            /*
            std::string debugFile = "vacInd_" + std::to_string(cnt)+".bp";
            Fcalc->output_indicator_function((fileManager.result_path()/"debug"/debugFile).c_str());
            cnt ++;
            */
        }
        
        if (forceCalculation.size() > 0){
            std::string debugFile = "vacInd_" + std::to_string(cnt);
            forceCalculation[0]->output_indicator_function((fileManager.result_path()/"debug"/debugFile).c_str());
        }
        

        // initialize solver and rhs
        auto linSolver = mag_tools::LU_solver<T>(a, L, bc, A_delta);

        // setup time stepping
        auto timeStepping = mag_tools::time_stepping(scen.timeStepping);

        // Newton-Raphson parameters
        auto solverParams = mag_tools::solver_parameters(scen.solverInput);
        auto iterMonitor = mag_tools::utils::iteration_monitor();

        const double eps = solverParams.relTol;
        double res = 1/eps;
        const int maxIter = solverParams.maxIter;
        const double minRes = solverParams.minRes;

        double startRes = 1.0;
        double relaxation = 1.0;
        int iter = 0;

        // time stepping
        while (!timeStepping.is_finished()){
            
            if (rank == 0){
                std::cout   << "\nTime step " << timeStepping.get_step()  + 1 
                            << " of " << timeStepping.number_of_steps() << std::endl;
            }

            // reset Newton monitoring variables
            res = 1/eps;
            iter = 0;
            iterMonitor.next_step();

            // update source terms and residual
            for(auto& src: fieldSources){
                src.update_source(timeStepping.get_time());
            }    
            
            // calculate residual for convergence monitoring
            startRes =linSolver.calc_orig_rhs_norm(zeroL);
            linSolver.assemble_rhs();

        
            // switch from relative -> absolute residual if necessary
            if (startRes < minRes){
                    startRes = 1.0;
            }
            
            // Newton iteration
            while ((res > eps) && (iter < maxIter)){

                iter ++;

                // solve linear system
                linSolver.assemble_matrix();
                linSolver.solve();
                
                // update old solution
                VecAXPY(_A.vec(), relaxation, _A_delta.vec());
                VecGhostUpdateBegin(_A.vec(), INSERT_VALUES, SCATTER_FORWARD);
                VecGhostUpdateEnd(_A.vec(), INSERT_VALUES, SCATTER_FORWARD);
                
                // update field quantities and dependent variables
                curlEval->update_result();
                for (auto& material: *materialDefinition){
                    material.update_material_models();
                }

                // assemble rhs and calculate current residual
                linSolver.assemble_rhs();
                res = linSolver.get_rhs_norm() / startRes;

                iterMonitor.next_iteration(res);
                if (rank == 0){
                    std::cout << "Iteration " << iter << " with residual: " << res << " and start "<< startRes << std::endl;
                }

            }

            // update field evolution (for hysteresis models)
            for (auto& material: *materialDefinition){
                material.prepare_next_timestep();
            }

            // update auxiliary variables
            magFieldCalc->update_result();

            // output
            B = dgProj.project_quad_to_dg(Bquad, B);
            H = dgProj.project_quad_to_dg(Hquad, H);
            M_out = dgProj.project_quad_to_dg(M, M_out);
            M0_out = dgProj.project_quad_to_dg(M0, M0_out);
            j_out= dgProjSca.project_quad_to_dg(j, j_out);

            
            for(auto& file: outputFiles){
                file->write(timeStepping.get_time());
            }

            // update force calculation            
            wmag->perform_update();
            for (auto& Fcalc:forceCalculation){
                Fcalc->update_forces(0.0);
                Fcalc->store_results(timeStepping.get_time());
            }
            
            // next time step
            timeStepping.next_time_step();   
        }


        // write additional output
        
        const std::string xmlFileName = (fileManager.result_path()/"results.xml").c_str();
        mag_tools::output::result_xml xmlOutput(xmlFileName);

        auto timeOutput = xmlOutput.append_time_series("time",{"unit"}, {"s"});
        timeOutput->add_values(timeStepping.get_time_vec());
        
        for (auto& Fcalc:forceCalculation){
           Fcalc->append_forces_to_output(xmlOutput);
        }

        xmlOutput.write_file();
        iterMonitor.write_statistic_file((fileManager.result_path()/"statistics.xml").c_str());
        
        if (rank == 0){    
            std::cout << "\n\nEnd of field calculation\n\n";
        }

    }
}