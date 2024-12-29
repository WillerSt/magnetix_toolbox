#include <magnetics_toolbox/coupled_material_model.h>

namespace mag_tools{

    coupled_atan_H_B_2D::coupled_atan_H_B_2D(
			const parameters_atan_H_B& parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr)
	:magnetic_material_coupling(VecDofsIn, JacDofsOut, VecDofsOut, blockNr),
			Bsat(parameters.Bsat), mumax(parameters.mumax), mulin(parameters.mulin)
	{
		Jacobean_H_B << 1/((mumax+mulin)*constants::mu0), 0, 0, 1/((mumax+mulin)*constants::mu0);
		this->calculate_jacobean();
	}


	// Update for scalar potential
	void coupled_atan_H_B_2D::update_coupling_entries_J_mu()
	{
		EigenVec J;


		H = {*(InputEntries[0][0]), *(InputEntries[1][0])};
		double Hnorm = H.norm();
		double Bnorm = H_B_curve(Hnorm);
		double Jnorm = Bnorm-constants::mu0*Hnorm;



		if (Hnorm > 1e-10)
		{
			B = Bnorm * H/Hnorm;
		}
		else
		{
			B = {0,0};
		}


		if (Hnorm > 1e-10)
		{
			J = Jnorm * H/Hnorm;
		}
		else
		{
			J = {0,0};
		}

		double jac11 = calc_Bx_Hx(H, Hnorm);
		double jac12 = calc_Bx_Hy(H, Hnorm);
		double jac21 = calc_By_Hx(H, Hnorm);
		double jac22 = calc_By_Hy(H, Hnorm);

		EigenJacob temp;
		temp << jac11, jac12, jac21, jac22; // This is necessary, find out why...
		Jacobean_H_B = temp;

		*(OutputMag[0][0])=J[0];
		*(OutputMag[1][0])=J[1];

		this->calculate_jacobean();

	}

	// update for vector potential
	void coupled_atan_H_B_2D::update_coupling_entries_M_nu()
	{
		EigenVec M;

		B = {*(InputEntries[0][0]), *(InputEntries[1][0])};
		double Bnorm = B.norm();
		double Hnorm = calc_H(Bnorm);
		double Mnorm = Bnorm*constants::mu0inv - Hnorm;

		if (Bnorm > 1e-10)
		{
			H = Hnorm * B/Bnorm;
		}
		else
		{
			H = {0,0};
		}


		if (Bnorm > 1e-10)
		{
			M = Mnorm * B/Bnorm;
		}
		else
		{
			M = {0,0};
		}

		double jac11 = calc_Bx_Hx(H, Hnorm);
		double jac12 = calc_Bx_Hy(H, Hnorm);
		double jac21 = calc_By_Hx(H, Hnorm);
		double jac22 = calc_By_Hy(H, Hnorm);

		EigenJacob temp;
		temp << jac11, jac12, jac21, jac22; // This is necessary, find out why...
		Jacobean_H_B = temp.inverse();

		*(OutputMag[0][0])=M[0];
		*(OutputMag[1][0])=M[1];
		this->calculate_jacobean();
	}
	/*
	std::vector<dolfin::la_index> get_associated_indices()
	{
		return VecDofs;
	}
	*/

	void coupled_atan_H_B_2D::set_jacobian_mu0()
	{
		*(OutputJac[0][0]) = constants::mu0;
		*(OutputJac[0][1])  = 0;
		*(OutputJac[1][0])  = 0;
		*(OutputJac[1][1])  = constants::mu0;
	}

	void coupled_atan_H_B_2D::reset_J(){
		EigenVec J = this->get_J();
		*(OutputMag[0][0])=J[0];
		*(OutputMag[1][0])=J[1];
	}

	void coupled_atan_H_B_2D::reset_M(){
		EigenVec M = this->get_M();
		*(OutputMag[0][0])=M[0];
		*(OutputMag[1][0])=M[1];
	}

	
	void coupled_atan_H_B_2D::calculate_jacobean()
	{
		*(OutputJac[0][0]) = Jacobean_H_B.coeff(0,0);
		*(OutputJac[0][1]) = Jacobean_H_B.coeff(0,1);
		*(OutputJac[1][0]) = Jacobean_H_B.coeff(1,0);
		*(OutputJac[1][1]) = Jacobean_H_B.coeff(1,1);
	}

	double coupled_atan_H_B_2D::H_B_curve(double Hnorm)
	{
		double Bnorm = constants::mu0 * (1+mulin) * Hnorm + fak1 * std::atan( fak2 * Hnorm);
		return Bnorm;
	}

	double coupled_atan_H_B_2D::calc_mur(double Hnorm)
	{
		return constants::mu0*(1+mulin) + constants::mu0 * (this->mumax -1) / (1+std::pow(fak2 * Hnorm,2));

	}

	double coupled_atan_H_B_2D::coupled_atan_H_B_2D::calc_Bx_Hx(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			return constants::mu0*(1+mulin) + fak1 * fak2  / (1+std::pow( (fak2*Hnorm), 2) ) * std::pow((Hx/Hnorm),2)
					+ fak1 * std::atan(fak2 * Hnorm) * std::pow((Hy/Hnorm),2) / Hnorm ;
		}
		else
		{
			return (this->mumax+mulin)*constants::mu0;
		}
	}

	double coupled_atan_H_B_2D::calc_By_Hy(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			return constants::mu0*(1+mulin) + fak1 * fak2 / (1+std::pow( (fak2*Hnorm), 2)) * std::pow((Hy/Hnorm),2)
					+ fak1  * std::atan(fak2 * Hnorm) * std::pow((Hx/Hnorm),2) / Hnorm ;
		}
		else
		{
			return (this->mumax+mulin)*constants::mu0;
		}
	}

	double coupled_atan_H_B_2D::calc_Bx_Hy(EigenVec H, double Hnorm)
	{
		double Hx = H[0];
		double Hy = H[1];
		if (Hnorm > 1e-6)
		{
			double fak3 = Hx * Hy / (std::pow(Hnorm,2));
			return fak1 * fak2 / (1+ std::pow((fak2 * Hnorm),2)) * fak3 - fak1 * std::atan(fak2*Hnorm) * fak3 / Hnorm;
		}
		else
		{
			return 0;
		}
	}
	double coupled_atan_H_B_2D::calc_By_Hx(EigenVec H, double Hnorm)
	{
	 return calc_Bx_Hy(H, Hnorm);
	}


	double coupled_atan_H_B_2D::calc_H (double Bnorm)
	{
		double mu, deltaH;
		double H = 0;
		double B = 0;
		double deltaB = Bnorm;
		double diff = 1;

		if (Bnorm > 1e-10)
		{
			while (diff > 1e-15)
			{
				mu = calc_mur(H);
				deltaH = deltaB/mu;
				H += deltaH;
				B = H_B_curve(H);
				deltaB = Bnorm-B;
				diff = std::abs(deltaB)/Bnorm;
			}
		}
		else
		{
			H=Bnorm/(constants::mu0*(mumax+mulin));
		}
		return H;
	}



}