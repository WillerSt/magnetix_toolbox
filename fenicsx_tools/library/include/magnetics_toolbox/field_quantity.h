#pragma once
#define USE_QUICK_VERSION

#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>
#include <curl_evaluation.h>
#include <magnetics_toolbox/mag_tools_basic.h>

#ifdef USE_QUICK_VERSION
#include <curl_projection.h>
#else
#include <linear_solver.h>
#endif

namespace mag_tools{
    template <typename T> class field_quantity{

        protected:
            std::shared_ptr<dolfinxFunction<T>> outFct;

        public:
        field_quantity(const std::shared_ptr<const dolfinxFS<T>>& funcSpaceIn)
            :outFct(std::make_shared<dolfinxFunction<T>>(dolfinxFunction<T>(funcSpaceIn))){};

        field_quantity(const std::shared_ptr<dolfinxFunction<T>>& funcIn):outFct(funcIn){};
        virtual ~field_quantity(){};

        virtual void update_result() = 0;

        std::shared_ptr<dolfinxFunction<T>> get_result_function(){
            return this->outFct;
        }
        
    };

    template <typename T> class linear_superposition: public field_quantity<T>{

        private:
        
        const std::vector<std::shared_ptr<dolfinxFunction<T>>> inFcts;
        const std::vector<double> scale;
        std::vector<dolfinxVector> dofVecs;

        dolfinxVector outVec;


        public:
        linear_superposition(const std::vector<std::shared_ptr<dolfinxFunction<T>>>& inputFunctions,  std::vector<double> scaleIn);
        void update_result();
    };

    template <typename T> class surf_curl_evaluation: public field_quantity<T>{


        private:
            const std::shared_ptr<const dolfinxFS<T>> quadFS;
            const std::shared_ptr<dolfinxFunction<T>> sol; 
            const T scale;

            const std::shared_ptr<const dolfinxConstant<T>> scaleConst = std::make_shared<const dolfinxConstant<T>>(scale);

            #ifdef USE_QUICK_VERSION
            const std::shared_ptr<varForm<T>> projForm = std::make_shared<varForm<T>>(dolfinx::fem::create_form<T>(
                *form_curl_projection_a, {quadFS, sol->function_space()}, {}, {{"c", scaleConst}}, {}, {}));
                
            dolfinxMatrix curlProjMat = dolfinx::la::petsc::Matrix(dolfinx::fem::petsc::create_matrix(*projForm), false);
            #else                    
                std::shared_ptr<varForm> evalForm = std::make_shared<varForm>(dolfinx::fem::create_form<T>(
                        *form_curl_evaluation_a, {quadFS, quadFS}, {std::vector<std::shared_ptr<const dolfinxFunction>>()}, {}, {}, {}));
                std::shared_ptr<varForm> evalRight = std::make_shared<varForm>(dolfinx::fem::create_form<T>(
                        *form_curl_evaluation_L, {quadFS, quadFS}, {std::vector<std::shared_ptr<const dolfinxFunction>>()}, {}, {}, {}));
                LU_solver<T> solver = LU_solver<T>(evalForm, evalRight, {}, this->outFct)
            #endif
            
            dolfinxVector curlQuadVec = dolfinxVector(dolfinx::la::petsc::create_vector_wrap(*this->outFct->x()), false);

            dolfinxVector solVec = dolfinxVector(dolfinx::la::petsc::create_vector_wrap(*sol->x()), false);

            dolfinxFunction<T> quadDiag = dolfinxFunction<T>(quadFS);

            dolfinxVector  quadDiagVec = dolfinxVector(dolfinx::la::petsc::create_vector_wrap(*quadDiag.x()), false);


        public:
        surf_curl_evaluation(const std::shared_ptr<dolfinxFunction<T>>& solIn, const std::shared_ptr<const dolfinxFS<T>>& quadFSIn, const T& scaleIn = 1.0);

        surf_curl_evaluation(const std::shared_ptr<dolfinxFunction<T>>& solIn, const std::shared_ptr<dolfinxFunction<T>>& quadFuncIn, const T& scaleIn = 1.0);
        void update_result();

        private:
        void initialize_variables();

        void update_quad_point_values();
    };
}

