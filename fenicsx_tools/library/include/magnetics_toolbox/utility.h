#pragma once
#include <magnetics_toolbox/mag_tools_basic.h>
namespace mag_tools::utility{
    template <typename T>
    std::shared_ptr<const dolfinxFS<T>> create_quadrature_functionspace(const std::shared_ptr<dolfinxMesh<T>>& mesh, const std::string& type){
        
        const int degQ = 2;
        const size_t dimG = mesh->topology()->dim();
        
        auto quadRule = basix::quadrature::make_quadrature<T>(basix::quadrature::type::Default, basix::cell::type::triangle, basix::polyset::type::standard, degQ);

        if (!strcmp(type.c_str(),"mat")){
            //std::cout << "Here mat\n";
            auto qMatE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, dimG*dimG));
            auto const dofmapMat = std::make_shared<const dolfinx::fem::DofMap>(
                            dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                            dolfinx::fem::create_element_dof_layout<T>(*qMatE), 
                            *mesh->topology(), nullptr, nullptr));

            auto quadMatFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qMatE, dofmapMat,{dimG,dimG}));
            return quadMatFS;
        }
        else if(!strcmp(type.c_str(),"vec")){
            //std::cout << "Here vec\n";
            auto qVecE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, dimG));
            
            auto const dofmapVec = std::make_shared<const dolfinx::fem::DofMap>(
                                dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                                dolfinx::fem::create_element_dof_layout<T>(*qVecE), 
                                *mesh->topology(), nullptr, nullptr));

            auto quadVecFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qVecE, dofmapVec,{dimG}));
            return quadVecFS;


        }
        else if(!strcmp(type.c_str(),"sca")){
            //std::cout << "Here sca\n";
            auto qScaE = std::make_shared<const dolfinx::fem::FiniteElement<T>>(dolfinx::fem::FiniteElement<T>(dolfinx::mesh::CellType::triangle, quadRule[0], {quadRule[0].size()/dimG, dimG}, 1));
            auto const dofmapSca = std::make_shared<const dolfinx::fem::DofMap>(
                    dolfinx::fem::create_dofmap(MPI_COMM_WORLD, 
                    dolfinx::fem::create_element_dof_layout<T>(*qScaE), 
                    *mesh->topology(), nullptr, nullptr));
            auto quadScaFS = std::make_shared<const dolfinx::fem::FunctionSpace<T>>(
                dolfinx::fem::FunctionSpace<T>(mesh, qScaE, dofmapSca,{1}));
            return quadScaFS;

        }
        else{
            std::cout << "ERROR: unknown type " << type << "for quadrature element\n";
            return nullptr;
        }

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
            const auto dgVecFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {dimG,1}));
            return dgVecFS;
        }
        else if(!strcmp(type.c_str(), "sca")){
            const auto dgScaFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {}));
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
            const auto cgVecFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {dimG,1}));
            return cgVecFS;
        }
        else if(!strcmp(type.c_str(), "sca")){
            const auto cgScaFS = std::make_shared<fem::FunctionSpace<T>>(dolfinx::fem::create_functionspace(mesh, dg1Elemem, {}));
            return cgScaFS;
        }
        else{
            std::cout << "ERROR: unknown type " << type << "for dg element\n";
            return nullptr;
        }
    }
}