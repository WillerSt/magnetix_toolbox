    #pragma once
    #include "magnetics_toolbox/mag_tools_basic.h"
    #include "magnetics_toolbox/scenario_description.h"

    namespace mag_tools{

        struct mesh_statistics{
            size_t nEdges;
            size_t nVertices;
            size_t nFaces;
            size_t nVolumes;
        };

        template<typename T> class mesh_container{

            private:
            const dolfinx::fem::CoordinateElement<U<T>> coordElement; 

            dolfinx::io::XDMFFile meshInput;
            dolfinx::io::XDMFFile boundaryInput;

            const std::shared_ptr<dolfinxMesh<T>> mesh;
            const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> meshMarkers;
            const std::shared_ptr<const dolfinx::mesh::MeshTags<std::int32_t>> boundaryMarkers;

            const int dimG = mesh->topology()->dim();

            public:
            mesh_container(const mag_tools::scenario_description& desc);

            public:
            auto inline get_mesh() const{
                return this->mesh;
            }
            auto inline get_mesh_markers() const{
                return this->meshMarkers;
            }
            auto inline get_boundary_markers() const{
                return this->boundaryMarkers;
            }
            private:
            std::shared_ptr<dolfinxMesh<T>> initialize_mesh();

        };
    }