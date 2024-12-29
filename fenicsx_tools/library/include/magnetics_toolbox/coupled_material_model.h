/*
 * magnetic_material_coupling.h
 *
 *  Created on: 22.09.2017
 *      Author: Stephan Willerich
 */

#pragma once
#include <dolfinx.h>
#include <magnetics_toolbox/quadrature_dof_coupler.h>
#include <Hysteresis_Operator/continuous_dpc_model.h>
#include <Hysteresis_Operator/xml_grid_constructor.h>

//#include "cxx_misc_tools/constants.h"
#include <math.h>

#include "magnetics_toolbox/constants.h"

namespace mag_tools{
// Interface class to hold parameters of arbitrary material model,
// every custom parameter class should inherit publicly form the interface.
class material_model_parameters
{
public:
	material_model_parameters(){};
	inline virtual ~material_model_parameters(){}
};

class akima_spline{
	const std::shared_ptr<const std::vector<double>> xTable;
	const std::shared_ptr<const std::vector<double>> yTable;

	const std::shared_ptr<const std::vector<double>> m;
	const std::shared_ptr<const std::vector<double>> s;

	const std::shared_ptr<const std::vector<double>> a;
	const std::shared_ptr<const std::vector<double>> b;
	const std::shared_ptr<const std::vector<double>> c;
	const std::shared_ptr<const std::vector<double>> d;

	const std::shared_ptr<const std::vector<double>> intYTable;

	const double xMax = (*xTable)[xTable->size()-1]; 
	const double xMin = (*xTable)[0]; 
	const double yStart = (*yTable)[0];
	const double yEnd = (*yTable)[xTable->size()-1];
	
	public:
	akima_spline(
		const std::shared_ptr<const std::vector<double>>& xTableIn, 
		const std::shared_ptr<const std::vector<double>>& yTableIn):
		xTable(xTableIn), yTable(yTableIn), 
		m(std::make_shared<const std::vector<double>>(calculate_m())),
		s(std::make_shared<const std::vector<double>>(calculate_s())),
		a(std::make_shared<const std::vector<double>>(calculate_a())),
		b(std::make_shared<const std::vector<double>>(calculate_b())),
		c(std::make_shared<const std::vector<double>>(calculate_c())),
		d(std::make_shared<const std::vector<double>>(calculate_d())),
		intYTable(std::make_shared<const std::vector<double>>(calculate_integral_table())){

		if (!perform_check(false)){
			std::cout << "WARNING: Makima inerpolation check failed\n";
		};

	}

	double interpolate_value(const double& x) const;

	std::pair<double, double> get_start_values(const double& yRef) const;
	
	double interpolate_derivative(const double& x) const;
	
	double interpolate_integral(const double& x) const;

	std::pair<double, double> interpolate_val_der(const double&  x) const;

	private:

	bool perform_check(bool displayValue);

	std::vector<double> calculate_m();
	std::vector<double> calculate_s();

	std::vector<double> calculate_a();
	std::vector<double> calculate_b();
	std::vector<double> calculate_c();
	std::vector<double> calculate_d();

	std::vector<double> calculate_integral_table();


};

//template <typename materialModel>
template <typename T> class magnetic_material_coupling
{ 

protected:
	//materialModel materialCharacteristic;

	// Coefficient indices in (global) FEA vectors
    using refVec = mag_tools::quadrature_dof_coupler<T>::refVecType;
    using refBlock  = std::vector<std::vector<T*>>;

	// Actual input and Output Values
	const refBlock InputEntries;
	const refBlock OutputMag;
	const refBlock OutputJac;


	const std::size_t blockNr;

	// Coefficient indices in coupling vectors (local, above)
	const bool JacobeanPresent = true;

	double energyDensity = 0.0;

	refBlock refEnergy; 

public:
	// Constructor with coupling vector and jacobean
	magnetic_material_coupling(
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNrIn)
			: InputEntries((*VecDofsIn)[blockNrIn]), OutputMag((*VecDofsOut)[blockNrIn]), OutputJac((*JacDofsOut)[blockNrIn]),
				blockNr(blockNrIn) 
	{};

	// constructor with coupling vector and no jacobean
	inline virtual ~magnetic_material_coupling(){}

protected:

	// Standard function to initialize Material parameters
	// virtual materialModel initializeMaterialModel(material_model_parameters parameters) = 0;

public:
	 // Interface functions
	 virtual void update_coupling_entries_M_nu() = 0;
	 virtual void update_coupling_entries_J_mu() = 0;
	 virtual void accept_state() = 0;
	 virtual void reset_M() = 0;
	 virtual void reset_J() = 0;

	 virtual void couple_energy_density(const std::shared_ptr<const refVec>& energyDensityDofsOut){
		refEnergy = (*energyDensityDofsOut)[blockNr];
	 }

	virtual void update_energy_density(){
		*(refEnergy[0][0]) =energyDensity;
	}
	 //virtual void gather_statistics(material_model_statistics& statistics) = 0;

};

class parameters_lin_H_B: public material_model_parameters{
	
	public:
	const double muR;

	public:
	parameters_lin_H_B(const double& muRIn):muR(muRIn){};



};

class parameters_lin_PM_2D: public material_model_parameters{
	
	
	public:
	using EigenVec =  continuous_dpc_model<2>::EigenType;
	const double muR;
	const double Br;
	const EigenVec dir;
	
	public:
	parameters_lin_PM_2D(const double& muRIn, const double BrIn, const double& xDir, const double& yDir):
	muR(muRIn), Br(BrIn), dir(initialize_direction(xDir, yDir)){};

	private:
	EigenVec initialize_direction(const double& xDir, const double& yDir){
		EigenVec direction = {xDir, yDir};
		direction.normalize();
		return direction;
	}



};

class parameters_atan_H_B 		: public material_model_parameters{

public:
	const double Bsat;
	const double mumax;
	const double mulin;


public:
	parameters_atan_H_B(double BsatIn, double mumaxIn)
		:Bsat(BsatIn),mumax(mumaxIn), mulin(0.0)
	{};
	parameters_atan_H_B(double BsatIn, double mumaxIn, double mulinIn)
		:Bsat(BsatIn),mumax(mumaxIn), mulin(mulinIn)
	{};
};

class parameters_DPC_Model_Cont_2D	: public material_model_parameters
{
public:
	const std::shared_ptr<const malloc_grid<2, CS_sphere_cont<2>>> dpc_grid;
	const continuous_dpc_model<2> defaultModel;

public:
	parameters_DPC_Model_Cont_2D(const std::shared_ptr<const malloc_grid<2, CS_sphere_cont<2>>>& dpc_gridIn)
			:dpc_grid(dpc_gridIn), defaultModel(dpc_grid)
	{
		// std::cout << "default model created\n";
	}
};

class parameters_spline_H_B: public material_model_parameters{	
	public:
	const std::shared_ptr<const std::vector<double>> Htable;
	const std::shared_ptr<const std::vector<double>> Btable;
	const std::shared_ptr<const akima_spline> interpolationFct;

	parameters_spline_H_B(const std::vector<double>& H, const std::vector<double>& B):
		Htable(std::make_shared<const std::vector<double>>(H)),
		Btable(std::make_shared<const std::vector<double>>(B)),
		interpolationFct(std::make_shared<const akima_spline>(akima_spline(Htable, Btable))){

	}
};

class coupled_atan_H_B_2D 		: public magnetic_material_coupling<PetscScalar>
{

public:
	using parameterType = parameters_atan_H_B;
	using EigenVec = Eigen::Matrix<double, 2, 1>;
	using EigenJacob = Eigen::Matrix<double, 2, 2>;

private:

	const double Bsat;
	const double mumax;
	const double mulin;

	const double fak1 = 2 / constants::pi * this->Bsat;
	const double fak2 = (constants::pi * constants::mu0 *(this->mumax -1)) / (2* this->Bsat);

	EigenVec B, H;

	EigenJacob Jacobean_H_B;

public:

	coupled_atan_H_B_2D(
			const parameters_atan_H_B& parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr);
	
	// Update for scalar potential
	void update_coupling_entries_J_mu();
	// update for vector potential
	void update_coupling_entries_M_nu();

	inline void accept_state(){};

	inline EigenVec get_J() const
	{
		return B - constants::mu0 * H;
	}

	inline EigenVec get_M() const
	{
		return B*constants::mu0inv - H;
	}

	inline EigenVec get_B() const
	{
		return B;
	}

	inline EigenVec get_H() const
	{
		return H;
	}

	/*
	std::vector<dolfin::la_index> get_associated_indices()
	{
		return VecDofs;
	}
	*/

	void set_jacobian_mu0();

	void reset_J();

	void reset_M();

	/*
	void gather_statistics(material_model_statistics& statistics)
	{
	}
	*/

private:
	void calculate_jacobean();

	double H_B_curve(double Hnorm);
	double calc_mur(double Hnorm);
	double calc_Bx_Hx(EigenVec H, double Hnorm);

	double calc_By_Hy(EigenVec H, double Hnorm);
	double calc_Bx_Hy(EigenVec H, double Hnorm);
	double calc_By_Hx(EigenVec H, double Hnorm);
	double calc_H (double Bnorm);
	

};

class coupled_DPC_Model_Cont_2D  : public magnetic_material_coupling<PetscScalar>
{

public:
	using parameterType = parameters_DPC_Model_Cont_2D;

private:
	typedef continuous_dpc_model<2>::EigenType EigenVec;
	typedef continuous_dpc_model<2>::EigenJacob EigenJacob;

	continuous_dpc_model<2> hystOp;

public:
	coupled_DPC_Model_Cont_2D(
			const parameters_DPC_Model_Cont_2D& parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr);



public:
	void update_coupling_entries_M_nu();

	void update_coupling_entries_J_mu();

	inline EigenVec get_J() const
	{
		return hystOp.get_J();
	}

	inline EigenVec get_M() const
	{
		return hystOp.get_M();
	}

	inline EigenVec get_B() const
	{
		return hystOp.get_B();
	}

	inline EigenVec get_H() const
	{
		return hystOp.get_H();
	}

	inline void accept_state()
	{
		hystOp.accept_evolution();
	}

	void display_info()
	{
		/*
		std::cout << "DPC Model coupled to \n" << "MagDofs: " << VecDofs[0] << " " << VecDofs[1]
			<< "\n" << "JacDofs: " << JacDofs[0] << " " << JacDofs[1] << " " << JacDofs[2] << " " << JacDofs[3] << "\n"
			<< "Current Input: " << (*InputEntries)[VecDofs[0]] << " " << (*InputEntries)[VecDofs[1]] << "\n"
			<< "Current Output: " << (*OutputMag)[VecDofs[0]] << " " << (*OutputMag)[VecDofs[1]] << "\n"
			<< "Current Output: " << hystOp.get_B()[0] << " " << hystOp.get_B()[1] << "\n";
			*/
	}


private:

	void calculate_jacobian_H_B();

	void calculate_jacobian_B_H();

    void calculate_jacobian_B_H_force_symmetry();

    void calculate_jacobian_H_B_force_symmetry();
public:

	void set_jacobian_mu0();

	void reset_J();

	void reset_M();
	/*
	void gather_statistics(material_model_statistics& statistics)
	{
		
		statistics.inversionIterations +=  hystOp.get_iteration_number();
		statistics.lineSearches +=  hystOp.get_line_search_number();
		statistics.lineSearchIterations +=  hystOp.get_line_search_iteration_number();
		statistics.failedInversions +=  hystOp.get_failed_inversions_number();
		statistics.failedLineSearches +=  hystOp.get_failed_line_search_number();
		statistics.maxDeviation = std::max(statistics.maxDeviation, hystOp.get_max_deviation());
		
	}*/
};

class coupled_lin_H_B : public magnetic_material_coupling<PetscScalar>{
	
	const double muR;
	const double mu;
	const double nu;
	const double mFac;
	
	public:
	coupled_lin_H_B(
			const parameters_lin_H_B parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr)
	:magnetic_material_coupling(VecDofsIn, JacDofsOut, VecDofsOut, blockNr), muR(parameters.muR), 
			mu(constants::mu0*this->muR), nu(1/mu), mFac((1- 1/muR)*(constants::mu0inv)){

		*(OutputJac[0][0]) = this->nu;
        *(OutputJac[0][1]) = 0.0;
        *(OutputJac[1][0]) = 0.0;
        *(OutputJac[1][1]) = this->nu;
	};

	public:
	 // Interface functions
	 void update_coupling_entries_M_nu(){
		*(OutputMag[0][0]) = *(InputEntries[0][0])*this->mFac;
		*(OutputMag[1][0]) = *(InputEntries[1][0])*this->mFac;		
	 };
	 void update_coupling_entries_J_mu(){};
	 void accept_state(){
		this->energyDensity = 0.5*this->nu*(std::pow(*(InputEntries[0][0]),2)+std::pow(*(InputEntries[1][0]),2));
	 };
	 void reset_M(){
		*(OutputMag[0][0]) = 0.0;
		*(OutputMag[1][0]) = 0.0;
	 };
	 void reset_J(){};
};

class coupled_const_mag : public magnetic_material_coupling<PetscScalar>{
	// TODO: Less wasteful implementation -> low priority
	const double nu;
	const double Br;
	const parameters_lin_PM_2D::EigenVec dir;
	
	public:
	coupled_const_mag(
			const parameters_lin_PM_2D parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr)	:magnetic_material_coupling(VecDofsIn, JacDofsOut, VecDofsOut, blockNr), 
		nu(1/(constants::mu0*parameters.muR)),Br(parameters.Br), dir(parameters.dir) {		
		*(OutputMag[0][0]) = dir[0]*Br*nu;
		*(OutputMag[1][0]) = dir[1]*Br*nu; 	
	};

	public:
	 // Interface functions
	 void update_coupling_entries_M_nu(){};
	 void update_coupling_entries_J_mu(){};
	 void accept_state(){};
	 void reset_M(){
		*(OutputMag[0][0]) = dir[0]*Br*nu;
		*(OutputMag[1][0]) = dir[1]*Br*nu;
	 };
	 void reset_J(){};
};

class coupled_spline_H_B_2D : public magnetic_material_coupling<PetscScalar>{
	using EigenVec = Eigen::Matrix<double, 2, 1>;
	using EigenMat = Eigen::Matrix<double, 2, 2>;

	private:
	const std::shared_ptr<const akima_spline> splineInt;

	EigenVec B = {0,0};
	EigenVec H = {0,0};

	double Hnorm = 0.0;

	double dHdB =  splineInt->interpolate_derivative(0.0);

	EigenMat Jacobean_H_B;
	
	public:
	coupled_spline_H_B_2D(
		const parameters_spline_H_B& parameters,
		const std::shared_ptr<const refVec>& VecDofsIn,
		const std::shared_ptr<const refVec>& JacDofsOut,
		const std::shared_ptr<const refVec>& VecDofsOut,
		const std::size_t& blockNr
	):magnetic_material_coupling(VecDofsIn, JacDofsOut, VecDofsOut, blockNr), splineInt(parameters.interpolationFct){
		Jacobean_H_B << splineInt->interpolate_derivative(0.0), 0.0, 0.0, splineInt->interpolate_derivative(0.0);
		this->calculate_jacobean();
	}
	public:
	// Interface functions
	 void update_coupling_entries_M_nu();

	void calculate_jacobean()
	{
		*(OutputJac[0][0]) = Jacobean_H_B.coeff(0,0);
		*(OutputJac[0][1]) = Jacobean_H_B.coeff(0,1);
		*(OutputJac[1][0]) = Jacobean_H_B.coeff(1,0);
		*(OutputJac[1][1]) = Jacobean_H_B.coeff(1,1);
	}
	 
	void update_coupling_entries_J_mu(){};
	 
	void accept_state(){
		this->energyDensity = H.norm()*B.norm()- splineInt->interpolate_integral(Hnorm);
		
		if (this->energyDensity < 0){
			std::cout << "Found it" << energyDensity << " H " << Hnorm << " B " << B.norm()  << " p " <<   H.norm()*B.norm() << " m "  << splineInt->interpolate_integral(Hnorm)<<std::endl;
		}
		
		//std::cout << "diff " <<  0.5*H.transpose()*B <<" "<<this->energyDensity << std::endl; 
	 };
	 
	void reset_M(){
		*(OutputMag[0][0]) = 0.0;
		*(OutputMag[1][0]) = 0.0;
	 };
	 
	void reset_J(){};

	 private:
	double calc_Bx_Hx(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			return splineInt->interpolate_derivative(abs(Hnorm)) * std::pow((Hx/Hnorm),2)
					+ splineInt->interpolate_value(abs(Hnorm))* std::pow((Hy/Hnorm),2) / Hnorm ;
		}
		else
		{
			return splineInt->interpolate_derivative(abs(Hnorm));
		}
	}

	double calc_By_Hy(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			return splineInt->interpolate_derivative(abs(Hnorm)) * std::pow((Hy/Hnorm),2)
					+ splineInt->interpolate_value(abs(Hnorm))* std::pow((Hx/Hnorm),2) / Hnorm ;
		}
		else
		{
			return splineInt->interpolate_derivative(abs(Hnorm));
		}
	}
	
	double calc_Bx_Hy(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			double fak3 = Hx * Hy / (std::pow(Hnorm,2));
			return splineInt->interpolate_derivative(abs(Hnorm)) * fak3 - splineInt->interpolate_value(abs(Hnorm)) * fak3 / Hnorm;
		}
		else
		{
			return 0;
		}
	}
	
	double calc_By_Hx(EigenVec H, double Hnorm)
	{
	 return calc_Bx_Hy(H, Hnorm);
	}

	double calc_H (double Bnorm);
};
}