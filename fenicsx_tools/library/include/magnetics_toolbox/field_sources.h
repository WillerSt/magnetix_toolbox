#pragma once
#include <magnetics_toolbox/quadrature_dof_coupler.h>

namespace mag_tools{

    namespace src{

    // interface for time variant (or constant) scalar function
    class time_variation{

        public:
        time_variation(){}

        virtual ~time_variation(){};

        virtual double calc_current(const double& t) const = 0;
        
    };

    // sinusoidal scalar function
    class sine_variation : public time_variation{

        private:
        const double amp;
        const double f;        
        const double phase;
        const double offset;

        public:
        sine_variation(const double& ampIn, const double& fIn, const double& phaseIn, const double& offsetIn);

        double calc_current(const double& t) const;

    };

    class value_table :  public time_variation{
        private: 
        const std::vector<double> time;
        const std::vector<double> values;
        
        std::size_t curStep = 0;
        const bool interpolate = false;

        public:
        value_table(const std::vector<double> timeIn, const std::vector<double>& valuesIn): time(timeIn), values(valuesIn){
            for(std::size_t i = 0; i<time.size()-1;i++){
                if (time[i]>=time[i+1]){
                    std::cout << "WARNING: time values are not monotonously increasing:" << "time["<<i<<"]: " << time[i]
                         << "time["<<i+1<<"]: " << time[i+1] << std::endl;
                }
            }
        };

        double calc_current(const double& t) const{
            double minDiff = 1e100;
            std::size_t minIndex = 0;
            double diff;
            double curVal;
            bool exact = false;
            
            
            for (std::size_t i = 0; i < time.size(); i++){
                //std::cout << "time in: " << t << "candidate: " << time[i] << std::endl;
                diff = std::fabs(time[i]-t);
                if (diff < 1e-10){
                    curVal = values[i];
                    exact = true;
                    break;              
                }
                if (diff < minDiff){
                    minDiff = diff;
                    minIndex = i;
                }
            }
            // linear interpolation
            if (exact==false){
                std::size_t i,j;
                if ((t-time[minIndex])>0){
                    i = minIndex;
                    j = minIndex+1;
                }
                else{
                    i = minIndex-1;
                    j = minIndex;
                }
                curVal = values[i] + (values[j]-values[i])/(t-time[i]);
            }
            return curVal;
        }


    };

    // class to model piecewise constant scalar source terms
    template<typename T> class const_scalar_source{

        using meshTags = dolfinx::mesh::MeshTags<std::int32_t>;

        private:
        std::shared_ptr<dolfinxFunction<T>> quadFunc;
        std::vector<quadrature_dof_coupler<T>> coupledDofs;
        const std::vector<double> scaleFactors;
        const std::vector<int> tagIdx; 

        const std::shared_ptr<const time_variation> curFunc;

        double curVal;


        
        public:
        const_scalar_source(const std::shared_ptr<dolfinxFunction<T>>& quadFuncIn, const std::shared_ptr<const meshTags>& meshMarkers , const std::vector<double>& scaleIn,
            const std::vector<int>& tagIdxIn, const std::shared_ptr<const time_variation>& curFuncIn);

        void update_source(const double& t);

        double get_excitation_value() const;

        void set_excitation_unity();

        void set_excitation_zero();

        private:
        void set_dofs_to_curVal();
    };

    }

}