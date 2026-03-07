// SPDX-FileCopyrightText: 2025 Stephan Willerich
// SPDX-License-Identifier: MIT License

#pragma once
#include "pugixml.hpp"
#include <magnetics_toolbox/mag_tools_basic.h>
namespace mag_tools::utility{


    template <typename T>
    concept DolfinxScalar = std::floating_point<T>;

    /*
    template <DolfinxScalar T>
    std::shared_ptr<const dolfinxFS<T>> create_quadrature_functionspace(const std::shared_ptr<dolfinxMesh<T>>& mesh, const std::string& type){
        
        const int degQ = 2;
        const size_t dimG = mesh->topology()->dim();
        
        auto quadRule = basix::quadrature::make_quadrature<T>(basix::quadrature::type::Default, basix::cell::type::triangle, basix::polyset::type::standard, degQ);

        if (!strcmp(type.c_str(),"mat")){
            //std::cout << "Here mat\n";
            auto qMatE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(
                dolfinx::mesh::CellType::triangle, 
                quadRule[0], 
                {quadRule[0].size()/dimG,dimG}, 
                {dimG, dimG}));
            auto const dofmapMat = std::make_shared<const dolfinx::fem::DofMap>(
                            dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                            dolfinx::fem::create_element_dof_layout<T>(*qMatE), 
                            *mesh->topology(), nullptr, nullptr));

            auto quadMatFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qMatE, dofmapMat));
            return quadMatFS;
        }
        else if(!strcmp(type.c_str(),"vec")){
            //std::cout << "Here vec\n";
            auto qVecE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(
                dolfinx::mesh::CellType::triangle, 
                quadRule[0], 
                {quadRule[0].size()/dimG,dimG}, 
                {dimG}));
            auto const dofmapVec = std::make_shared<const dolfinx::fem::DofMap>(
                                dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                                dolfinx::fem::create_element_dof_layout<T>(*qVecE), 
                                *mesh->topology(), nullptr, nullptr));

            auto quadVecFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qVecE, dofmapVec));
            return quadVecFS;


        }
        else if(!strcmp(type.c_str(),"sca")){
            //std::cout << "Here sca\n";
            auto qScaE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}));
            auto const dofmapSca = std::make_shared<const dolfinx::fem::DofMap>(
                    dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                    dolfinx::fem::create_element_dof_layout<T>(*qScaE), 
                    *mesh->topology(), nullptr, nullptr));
            auto quadScaFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qScaE, dofmapSca));
            return quadScaFS;

        }
        else{
            std::cout << "ERROR: unknown type " << type << "for quadrature element\n";
            return nullptr;
        }

    }
    */


    template <DolfinxScalar T>
    std::shared_ptr<const dolfinxFS<T>> create_quadrature_functionspace(const std::shared_ptr<dolfinxMesh<T>>& mesh, const std::string& type){
        const int degQ = 2;
        const size_t dimG = mesh->topology()->dim();
        auto basix_cell = dolfinx::mesh::cell_type_to_basix_type(mesh->topology()->cell_type());

        auto quadRule = basix::quadrature::make_quadrature<T>(basix::quadrature::type::Default,
                                                               basix_cell,
                                                               basix::polyset::type::standard,
                                                               degQ);
        auto &qpts = quadRule[0];
        const std::size_t n_q = qpts.size() / dimG;


        auto make_fs = [&](std::vector<std::size_t> value_shape){
            // pshape is the quadrature point shape: {n_q, dimG}
            std::array<std::size_t, 2> pshape = {n_q, dimG};

            auto fe = std::make_shared<const dolfinx::fem::FiniteElement<T>>(
                dolfinx::fem::FiniteElement<T>(mesh->topology()->cell_type(),
                                               qpts,
                                               pshape,
                                               value_shape));
            auto dofmap = std::make_shared<const dolfinx::fem::DofMap>(
                dolfinx::fem::create_dofmap(MPI_COMM_WORLD,
                                            dolfinx::fem::create_element_dof_layout<T>(*fe),
                                            *mesh->topology(), nullptr, nullptr));
            return std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, fe, dofmap));
        };

        if (type == "mat")
            return make_fs({dimG, dimG});
        if (type == "vec")
            return make_fs({dimG});
        if (type == "sca")
            return make_fs({});

        std::cerr << "ERROR: unknown type " << type << " for quadrature element\n";
        return nullptr;
    }
    
    template <typename T>
    std::shared_ptr<const dolfinxFS<T>> create_dg_functionspace(const std::shared_ptr<dolfinxMesh<T>>& mesh, const int& deg, const std::string& type = "sca"){
        
        const size_t dimG = mesh->topology()->dim();
        const auto cType = dolfinx::mesh::cell_type_to_basix_type(mesh->topology()->cell_type());

        auto dg1Elemem = basix::create_element<U<T>>(
            basix::element::family::P, cType, deg,
            basix::element::lagrange_variant::unset,
            basix::element::dpc_variant::unset, true);  

  

        if (!strcmp(type.c_str(), "vec")){
            const auto dgVecFS = std::make_shared<dolfinx::fem::FunctionSpace<U<T>>>(dolfinx::fem::create_functionspace<U<T>>(
                mesh, std::make_shared<const dolfinx::fem::FiniteElement<U<T>>>(dg1Elemem, std::vector<size_t>({dimG,1}))));
            return dgVecFS;
        }
        else if(!strcmp(type.c_str(), "sca")){
            const auto dgScaFS = std::make_shared<dolfinx::fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace<U<T>>(mesh, 
                std::make_shared<const dolfinx::fem::FiniteElement<U<T>>>(dg1Elemem)));  
            return dgScaFS;
        }
        else{
            std::cout << "ERROR: unknown type " << type << "for dg element\n";
            return nullptr;
        }
    }


    template <typename T>
    std::shared_ptr<const dolfinxFS<T>> create_cg_functionspace(const std::shared_ptr<dolfinxMesh<T>>& mesh, const int& deg, const std::string& type = "sca"){
        
        const size_t dimG = mesh->topology()->dim();
        const auto cType = dolfinx::mesh::cell_type_to_basix_type(mesh->topology()->cell_type());

        auto dg1Elemem = basix::create_element<U<T>>(
            basix::element::family::P, cType, deg,
            basix::element::lagrange_variant::unset,
            basix::element::dpc_variant::unset, false);  

        if (!strcmp(type.c_str(), "vec")){
            const auto cgVecFS = std::make_shared<dolfinx::fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace<U<T>>(mesh, 
                std::make_shared<const dolfinx::fem::FiniteElement<U<T>>>(dg1Elemem, {dimG,1})));
            return cgVecFS;
        }
        else if(!strcmp(type.c_str(), "sca")){
            auto cgScaFS = std::make_shared<dolfinx::fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace<U<T>>(mesh, 
                std::make_shared<const dolfinx::fem::FiniteElement<U<T>>>(dg1Elemem)));
            return cgScaFS;
        }
        else{
            std::cout << "ERROR: unknown type " << type << "for dg element\n";
            return nullptr;
        }
    }

    class iteration_monitor {
        std::vector<int> seriesSteps = {};
        std::vector<std::vector<double>> iterations = {};

        pugi::xml_document xmlDoc = pugi::xml_document();
        pugi::xml_node desc = xmlDoc.append_child(pugi::node_declaration);
        pugi::xml_node iterStatNode = xmlDoc.append_child("IterationStatistic");

        const int rank;

    public:
        iteration_monitor() : rank(get_rank()) {
            desc.append_attribute("version") = "1.0";
        }

        void next_step() {
            iterations.push_back({});
        }

        void next_iteration(const double& res) {
            iterations.back().push_back(res);
        }

        void write_statistic_file(const std::string& fileName) {
            if (rank == 0) {
                try {
                    for (std::size_t step = 0; step < iterations.size(); ++step) {
                        auto iterationsNode = iterStatNode.append_child("Iteration");
                        iterationsNode.append_attribute("step") = step;
                        
                        // Pre-allocate string to avoid repeated reallocations
                        std::string textContent;
                        textContent.reserve(iterations[step].size() * 16);  // rough estimate
                        
                        for (std::size_t i = 0; i < iterations[step].size(); ++i) {
                            textContent += " ";
                            textContent += std::to_string(iterations[step][i]);
                        }
                        iterationsNode.text().set(textContent.c_str());
                    }
                    
                    bool save_ok = xmlDoc.save_file(fileName.c_str());
                    if (!save_ok) {
                        std::cerr << "ERROR: Failed to write " << fileName << "\n";
                    }
                } catch (const std::exception& e) {
                    std::cerr << "ERROR in write_statistic_file: " << e.what() << "\n";
                }
            }
            
            // Only rank 0 writes; no need to barrier all ranks (they're not doing I/O)
            // If you must synchronize for safety, use MPI_Barrier after rank 0 finishes
            MPI_Barrier(MPI_COMM_WORLD);  // safe: ensures file is written before any rank continues
        }

    private:
        int get_rank() const {
            int currentRank;
            MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);
            return currentRank;
        }
    };
}