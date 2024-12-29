#include <magnetics_toolbox/coupled_material_model.h>

namespace mag_tools{   

    coupled_DPC_Model_Cont_2D::coupled_DPC_Model_Cont_2D(
			const parameters_DPC_Model_Cont_2D& parameters,
			const std::shared_ptr<const refVec>& VecDofsIn,
			const std::shared_ptr<const refVec>& JacDofsOut,
            const std::shared_ptr<const refVec>& VecDofsOut,
            const std::size_t& blockNr)
	:magnetic_material_coupling(VecDofsIn, JacDofsOut, VecDofsOut, blockNr), hystOp(parameters.defaultModel)
	{
		this->calculate_jacobian_H_B();
	}

    void coupled_DPC_Model_Cont_2D::update_coupling_entries_M_nu()
	{
		EigenVec NewIn;
		hystOp.calculate_evolution_B(EigenVec {*(InputEntries[0][0]), *(InputEntries[1][0])});
		NewIn = hystOp.get_M();
		*(OutputMag[0][0])=NewIn[0];
		*(OutputMag[1][0])=NewIn[1];
		this->calculate_jacobian_H_B();
	}

	void coupled_DPC_Model_Cont_2D::update_coupling_entries_J_mu()
	{
		EigenVec NewIn;
		hystOp.calculate_evolution_H(EigenVec {*(InputEntries[0][0]), *(InputEntries[1][0])});
		NewIn = hystOp.get_J();
		*(OutputMag[0][0])=NewIn[0];
		*(OutputMag[1][0])=NewIn[1];
		this->calculate_jacobian_B_H();
	}
    
    void coupled_DPC_Model_Cont_2D::calculate_jacobian_H_B()
    {
        EigenJacob LocalJacobian = hystOp.approximate_jacobian_H_B();
        *(OutputJac[0][0]) = LocalJacobian.coeff(0,0);
        *(OutputJac[0][1]) = LocalJacobian.coeff(0,1);
        *(OutputJac[1][0]) = LocalJacobian.coeff(1,0);
        *(OutputJac[1][1]) = LocalJacobian.coeff(1,1);
    }

    void coupled_DPC_Model_Cont_2D::calculate_jacobian_B_H()
    {
        EigenJacob LocalJacobian = hystOp.approximate_jacobian_B_H_explicitly();
        *(OutputJac[0][0]) = LocalJacobian.coeff(0,0);
        *(OutputJac[0][1]) = LocalJacobian.coeff(0,1);
        *(OutputJac[1][0]) = LocalJacobian.coeff(1,0);
        *(OutputJac[1][1]) = LocalJacobian.coeff(1,1);
    }

    void coupled_DPC_Model_Cont_2D::calculate_jacobian_B_H_force_symmetry()
    {
        EigenJacob LocalJacobian = hystOp.approximate_jacobian_B_H_explicitly();
        *(OutputJac[0][0]) = LocalJacobian.coeff(0,0);
        *(OutputJac[0][1]) =  0.5*(LocalJacobian.coeff(0,1)+LocalJacobian.coeff(1,0));
        *(OutputJac[1][0]) =  0.5*(LocalJacobian.coeff(0,1)+LocalJacobian.coeff(1,0));
        *(OutputJac[1][1]) = LocalJacobian.coeff(1,1);
    }

    void coupled_DPC_Model_Cont_2D::calculate_jacobian_H_B_force_symmetry()
    {
        EigenJacob LocalJacobian = hystOp.approximate_jacobian_H_B();
        *(OutputJac[0][0]) = LocalJacobian.coeff(0,0);
        *(OutputJac[0][1]) =  0.5*(LocalJacobian.coeff(0,1)+LocalJacobian.coeff(1,0));
        *(OutputJac[1][0]) =  0.5*(LocalJacobian.coeff(0,1)+LocalJacobian.coeff(1,0));
        *(OutputJac[1][1]) = LocalJacobian.coeff(1,1);
    }

    void coupled_DPC_Model_Cont_2D::set_jacobian_mu0()
    {
        *(OutputJac[0][0]) = constants::mu0;
        *(OutputJac[0][1])  = 0;
        *(OutputJac[1][0])  = 0;
        *(OutputJac[1][1])  = constants::mu0;
    }

    void coupled_DPC_Model_Cont_2D::reset_J(){
        EigenVec J = this->get_J();
        *(OutputMag[0][0])=J[0];
        *(OutputMag[1][0])=J[1];
    }

    void coupled_DPC_Model_Cont_2D::reset_M(){
        EigenVec M = this->get_M();
        *(OutputMag[0][0])=M[0];
        *(OutputMag[1][0])=M[1];
    }
}