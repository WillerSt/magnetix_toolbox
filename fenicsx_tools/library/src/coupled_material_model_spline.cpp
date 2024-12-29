#include "magnetics_toolbox/coupled_material_model.h"
namespace mag_tools{	
    double akima_spline::interpolate_value(const double& x) const{
		double y;
		if (x >= xMax){
			y = (yEnd + s->back() * (x-xMax));
		}
		else if(x<=xMin) {
			y = (yStart + (*s)[0] * (x-xMin));
		}
		else{
			for (std::size_t i = 0; i < xTable->size(); i++){
				if (i == xTable->size()-1){
					std::cout << "Should not get here akima interpolation x:"<< x << " xMin " << xMin << " xMax " << xMax <<  " \n";
					exit(1);
				}
				if (((*xTable)[i] <= x) && ((*xTable)[i+1]>x)){
					double dx = x - (*xTable)[i];
					y = (*a)[i] + (*b)[i]*dx + (*c)[i]*pow(dx,2) + (*d)[i]*pow(dx,3); 
					break;
				}

			}
		}
		return y;
	}

	std::pair<double, double> akima_spline::get_start_values(const double& yRef) const {
		std::size_t i;
			for (i = 0; i < yTable->size(); i++){
				if (((*yTable)[i] <= yRef) && ((*yTable)[i+1]>yRef)){
					break;
				}
			}
		
		return std::pair<double, double> ((*xTable)[i], (*yTable)[i]);

	}


	double akima_spline::interpolate_derivative(const double& x) const{
		double dy;
		if (x >= xMax){
			dy = s->back();
		}
		else if(x<=xMin) {
			dy = (*s)[0];
		}
		else{
			for (std::size_t i = 0; i < xTable->size(); i++){
				if (((*xTable)[i] <= x) && ((*xTable)[i+1]>x)){
					double dx = x - (*xTable)[i];
					dy = (*b)[i] + 2*(*c)[i]*dx + 3*(*d)[i]*pow(dx,2);
					if (dy < 0){
						std::cout << "alarm x: " << x  << " i: " << i << " dy: " << dy << "\n"; 
						std::cout << "alarm b: " << (*b)[i]  << " c: " << (*c)[i] << " d: " << (*d)[i] << "\n"; 
						dy = (*m)[i];
					}
					break;
				}
			}
		}
		return dy;
	}
	
	double akima_spline::interpolate_integral(const double& x) const{
		double yI;
		if (x > xMax){
			yI = intYTable->back() + yEnd * (x-xMax) + s->back() * (pow(x-xMax, 2));
		}
		else if(x<xMin) {
			yI = intYTable->front() + yStart * (xMin-x) + (*s)[0] *0.5* (pow(xMin-x,2));
		}
		else{
			for (std::size_t i = 0; i < xTable->size(); i++){
				if (i == xTable->size()-1){
					std::cout << "Should not get here akima interpolation x:"<< x << " xMin " << xMin << " xMax " << xMax <<  " \n";
					exit(1);
				}
				if (((*xTable)[i] <= x) && ((*xTable)[i+1]>x)){
					double xUp = x;
					double xLow = (*xTable)[i];
					yI = (*intYTable)[i] + (*a)[i]*(xUp-xLow) + 0.5*(*b)[i]*(pow(xUp-xLow,2)) 
						+ 1/3*(*c)[i]*(pow(xUp-xLow,3)) + 0.25*(*d)[i]*(pow(xUp-xLow,4));

					break;
				}

			}
		}
		return yI;
	}

	std::pair<double, double> akima_spline::interpolate_val_der(const double&  x) const {
		double y;
		double dy;
		if (x > xMax){
			y = (yEnd + s->back() * (x-xMax));
			dy = s->back();
		}
		else if(x<xMin) {
			y = (yStart + (*s)[0] * (x-xMin));
			dy = (*s)[0] ;
		}
		else{
			for (std::size_t i = 0; i < xTable->size(); i++){
				if (((*xTable)[i] <= x) && ((*xTable)[i+1]>x)){
					double dx = x - (*xTable)[i];
					y = (*a)[i] + (*b)[i]*dx + (*c)[i]*pow(dx,2) + (*d)[i]*pow(dx,3); 
					dy = (*b)[i] + 2*(*c)[i]*dx + 3*(*d)[i]*pow(dx,2);
					break;
				}
			}
		}
		return {y, dy};

	}

	bool akima_spline::perform_check(bool displayValue){
		bool checkPassed = true;
		for (std::size_t i = 0; i < xTable->size(); i++){
			auto yInt = interpolate_value((*xTable)[i]);
			checkPassed = checkPassed && abs((*yTable)[i]  - yInt) < 1e-6;
			if (displayValue == true){
				
				std::cout << "x: " << (*xTable)[i] << "; y: " << (*yTable)[i] << "; yInt: " << yInt 
						<< " diff " << (*yTable)[i]  - yInt << std::endl;
			}
			
		}
		return checkPassed;
	}

	std::vector<double> akima_spline::calculate_m(){
		std::vector<double> m;
		for (std::size_t i = 0; i < xTable->size()-1; i++){
			m.push_back( ((*yTable)[i+1]-(*yTable)[i])/ ((*xTable)[i+1]-(*xTable)[i]));
		}
		return m;
	}
	std::vector<double> akima_spline::calculate_s(){
		std::vector<double> s;
		for (std::size_t i = 0; i < xTable->size(); i++){
			if (i == 0){
				s.push_back((*m)[0]);
			}
			else if(i == 1){
				s.push_back(((*m)[1]+(*m)[0])*0.5);
			}
			else if(i == xTable->size()-2){
				s.push_back(((*m)[i-2]+(*m)[i-1])*0.5);
			}
			else if(i == xTable->size()-1){
				s.push_back((*m)[i-1]);
			}
			else{
				double den = abs((*m)[i+1]-(*m)[i])+fabs((*m)[i-1]-(*m)[i-2]);
				if (den < 1e-9){
					s.push_back(((*m)[i]+(*m)[i-1])*0.5);
				}
				else{
					s.push_back(
						(abs((*m)[i+1]-(*m)[i])*(*m)[i-1] + abs((*m)[i-1]-(*m)[i-2])*(*m)[i]) / den
					);
				}
				
			}
		}
		return s;
	}

	std::vector<double> akima_spline::calculate_a(){
		std::vector<double> a;
		for (size_t i = 0; i<yTable->size()-1; i++){
			a.push_back((*yTable)[i]);
		}
		return a;
	}

	std::vector<double> akima_spline::calculate_b(){
		std::vector<double> b;
		for (size_t i = 0; i<yTable->size()-1; i++){
			b.push_back((*s)[i]);
		}
		return b;
		
	}	
	
	std::vector<double> akima_spline::calculate_c(){
		std::vector<double> c;
		for (std::size_t i = 0; i<yTable->size()-1; i++){
			c.push_back((3*(*m)[i]-2*(*s)[i]-(*s)[i+1])/((*xTable)[i+1]-(*xTable)[i]));
		}
		return c;
	}

	std::vector<double> akima_spline::calculate_d(){
		std::vector<double> d;
		for (std::size_t i = 0; i<yTable->size()-1; i++){
			d.push_back( ((*s)[i]+(*s)[i+1]-2*(*m)[i]) / std::pow(((*xTable)[i+1]-(*xTable)[i]),2));
		}
		return d;		
	}
	
	std::vector<double> akima_spline::calculate_integral_table(){
		std::cout << "Calculating table\n";
		std::vector<double> yInt;
		yInt.push_back(0.0);
		double xUp, xLow;
		for (std::size_t i = 1; i<yTable->size(); i++){
			xUp = (*xTable)[i];
			xLow = (*xTable)[i-1];
			yInt.push_back(yInt[i-1] + (*a)[i]*(xUp-xLow) + 0.5*(*b)[i]*(pow(xUp-xLow,2)) 
				+ 1/3*(*c)[i]*(pow(xUp-xLow,3)) + 0.25*(*d)[i]*(pow(xUp-xLow,4)) );

			if (yInt[i] > (*xTable)[i] * (*yTable)[i])
			{
				std::cout << "Nope!: " << " H " << (*xTable)[i] << " B "<< (*yTable)[i]  << " yI " << yInt[i] << std::endl;
			}

			// std::cout << "xLow " << xLow << " xUp "  << xUp << " result: " << yInt[i]<< "\n";
			// std::cout << "a" << (*a)[i] << " b "  << (*b)[i] << " c " << (*c)[i]<< " d " <<(*d)[i]<< "\n";
		}
		return yInt;
	}

    void coupled_spline_H_B_2D::update_coupling_entries_M_nu(){
		EigenVec M;

		B = {*(InputEntries[0][0]), *(InputEntries[1][0])};
		double Bnorm = B.norm();
		Hnorm = calc_H(Bnorm);
		double Mnorm = Bnorm*constants::mu0inv - Hnorm;

		if (Bnorm > 1e-10)
		{
			H = Hnorm * B/Bnorm;
		}
		else
		{
			H = {0,0};
			Hnorm = 0.0;
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

		EigenMat temp;
		temp << jac11, jac12, jac21, jac22; // This is necessary, find out why...
		Jacobean_H_B = temp.inverse();

		*(OutputMag[0][0])=M[0];
		*(OutputMag[1][0])=M[1];
		this->calculate_jacobean();	
	 }

     	double coupled_spline_H_B_2D::calc_H (double Bnorm){		
		
		double mu, deltaH;

		auto startVals = splineInt->get_start_values(Bnorm);
		double Hnew = startVals.first;
		double Bnew = startVals.second;
		double deltaB = Bnorm-Bnew;
		double diff = std::abs(deltaB)/Bnorm;

		int iter = 0;
		double relaxation = 1.0;

		if (Bnorm > 1e-10)
		{
			while (diff > 1e-10)
			{	iter += 1;
				mu = splineInt->interpolate_derivative(abs(Hnew));
				deltaH = deltaB/mu;
				Hnew += deltaH*relaxation;
				Bnew = ((Hnew>0)-(Hnew<0))*splineInt->interpolate_value(abs(Hnew));
				deltaB = Bnorm-Bnew;
				diff = std::abs(deltaB)/Bnorm;
				
				if ( diff > 0.1){
					relaxation = 0.25;
				}else{
					relaxation = 1;
				}
					

				if (iter > 40){
					std::cout << "B: " << Bnew << " H: "<< Hnew  << " muR: "<<  mu/ constants::mu0 <<std::endl;
					if (iter > 50){
						std::cout << "Too much "<<diff<< "\n";
						exit(1);
					}
				}
			}
		}
		else
		{
			Hnew=Bnorm/(splineInt->interpolate_derivative(0.0));
		}
		//std::cout << " Success B: " << Bnorm << " H: " << Hnew  << " iter: " << iter <<std::endl;
		return Hnew;
	}

}