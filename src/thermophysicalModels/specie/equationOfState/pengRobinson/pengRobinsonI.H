/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#include "pengRobinson.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{





// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline scalar pengRobinson::pcrit()const
{
return pcrit_;
}



inline scalar pengRobinson::Tcrit()const
{
return Tcrit_;
}


inline scalar pengRobinson::rhostd()const
{
return rhostd_;
}


// Returns the Azentric Factor (Acentric Factor)
inline scalar pengRobinson::azentricFactor() const
{
return azentricFactor_;
}

//returns the pressure for a given density and temperature    
inline scalar pengRobinson::pReturn(const scalar rho,const scalar T) const
{
return  this->RR*T/((this->W()/rho)-this->b())-(this->a()*pow((1+this->n()*(1-pow((T/this->Tcrit()),0.5))),2))/(pow((this->W()/rho),2)+2*this->b()*(this->W()/rho)-pow(this->b(),2));
}



// Factor a of the pengRobinson equation of state  
//(molar values)
inline scalar pengRobinson::a() const
{
return 0.457235*pow(this->RR,2)*pow(this->Tcrit(),2)/this->pcrit();
}



// Factor b of the pengRobinson equation of state 
//(molar values)
inline scalar pengRobinson::b() const
{
return 0.077796*this->RR*this->Tcrit()/this->pcrit(); 
}

// Factor nb of the pengRobinson equation of state 
//(molar values)
inline scalar pengRobinson::n() const
{
// this equation is only valid for asentricFactors <0.49. 
return 0.37464+1.54226*this->azentricFactor()-0.26992*pow(this->azentricFactor(),2); 
}


//* * * * * * * * * * * * * Derivatives * * * * * * * * * * * //



//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar pengRobinson::dpdv(const scalar rho, const scalar T) const
{
return -(4*this->a()*this->n()*this->Tcrit()*(this->b()-(this->W()/rho))*(pow(this->b(),2)-pow((this->W()/rho),2))*(this->n()+1)
*pow((T/this->Tcrit()),0.5)
+this->Tcrit()*(-2*this->a()*pow((this->n()+1),2)*(pow(this->b(),3)-pow(this->b(),2)*(this->W()/rho)
-this->b()*pow((this->W()/rho),2)+pow((this->W()/rho),3))
+this->RR*T*(pow(this->b(),4)-4*pow(this->b(),3)*(this->W()/rho)+2*pow(this->b(),2)*pow((this->W()/rho),2)
+4*this->b()*pow((this->W()/rho),3)+pow((this->W()/rho),4)))
-2*this->a()*pow(this->n(),2)*T*(pow(this->b(),3)-pow(this->b(),2)*(this->W()/rho)-this->b()*pow((this->W()/rho),2)+pow((this->W()/rho),3)))
/(this->Tcrit()*pow((this->b()-(this->W()/rho)),2)*pow((pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)),2));
}



//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar pengRobinson::dpdT(const scalar rho, const scalar T) const
{
return (-this->a()*this->n()*(this->n()+1)*pow((T/this->Tcrit()),0.5)/(T*(pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)))+this->a()*pow(this->n(),2)/(this->Tcrit()*(pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)))-this->RR/(this->b()-(this->W()/rho)));
}




//Real  deviative dv/dT at constant pressure 
//by using implicit differentiation
//(molar values)
inline scalar pengRobinson::dvdT(const scalar rho,const scalar T) const
{
return (-1)*this->dpdT(rho,T)/this->dpdv(rho,T);
}




//(molar values)
inline scalar pengRobinson::dvdp(const scalar rho,const scalar T)  const 
{
return 1/this->dpdv(rho,T);
}





//(molar values)
//needed to calculate the internal energy
inline scalar pengRobinson::integral_p_dv(const scalar rho,const scalar T) const
{
return  pow(2,0.5)*this->a()*(2*this->n()*this->Tcrit()*(this->n()+1)*pow(T/this->Tcrit(),0.5)
-this->Tcrit()*(pow(this->n(),2)+2*this->n()+1)-pow(this->n(),2)*T)*log(this->b()*(1-pow(2,0.5))+(this->W()/rho))/(4*this->b()*this->Tcrit())
+this->RR*T*log((this->W()/rho)-this->b())
-pow(2,0.5)*this->a()*(2*this->n()*this->Tcrit()*(this->n()+1)*pow(T/this->Tcrit(),0.5)
-this->Tcrit()*(pow(this->n(),2)+2*this->n()+1)-pow(this->n(),2)*T)
*log(this->b()*(pow(2,0.5)+1)+(this->W()/rho))/(4*this->b()*this->Tcrit());
}


//(molar values)
//needed to calculate the entropy
inline scalar pengRobinson::integral_dpdT_dv(const scalar rho,const scalar T) const
{

return (pow(2,0.5)*this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)/(4*this->b()*T)
-pow(2,0.5)*this->a()*pow(this->n(),2)/(4*this->b()*this->Tcrit()))
*log(this->b()*(1-pow(2,0.5))+(this->W()/rho))
+this->RR*log((this->W()/rho)-this->b())
+(pow(2,0.5)*this->a()*pow(this->n(),2)/(4*this->b()*this->Tcrit())
-pow(2,0.5)*this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)/(4*this->b()*T))
*log(this->b()*(pow(2,0.5)+1)+(this->W()/rho));

}




//* * * * * * * * * * * * * second order Derivative based functions * * * * * * * * * * * //



//(molar values)
inline scalar pengRobinson::d2pdT2(const scalar rho,const scalar T) const
{
return this->a()*this->n()*(this->n()+1)*pow((T/this->Tcrit()),0.5)/(2*pow(T,2)*(pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)));
}




//(molar values)
inline scalar pengRobinson::d2pdv2(const scalar rho,const scalar T) const
{
return 2*(2*this->a()*this->n()*this->Tcrit()*(this->b()-(this->W()/rho))*(this->n()+1)*
(5*pow(this->b(),4)-4*pow(this->b(),3)*(this->W()/rho)-4*pow(this->b(),2)*pow((this->W()/rho),2)
+3*pow((this->W()/rho),4))*pow(T/this->Tcrit(),0.5)
+this->Tcrit()*(this->RR*T*(pow(this->b(),6)-6*pow(this->b(),5)*(this->W()/rho)+9*pow(this->b(),4)*pow((this->W()/rho),2)
+4*pow(this->b(),3)*pow((this->W()/rho),3)-9*pow(this->b(),2)*pow((this->W()/rho),4)-6*this->b()*pow((this->W()/rho),5)
-pow((this->W()/rho),6))
-this->a()*pow((this->n()+1),2)*(5*pow(this->b(),5)-9*pow(this->b(),4)*(this->W()/rho)
+4*pow(this->b(),2)*pow((this->W()/rho),3)+3*this->b()*pow((this->W()/rho),4)-3*pow((this->W()/rho),5)))
-this->a()*pow(this->n(),2)*T*(5*pow(this->b(),5)-9*pow(this->b(),4)*(this->W()/rho)+4*pow(this->b(),2)*pow((this->W()/rho),3)
+3*this->b()*pow((this->W()/rho),4)-3*pow((this->W()/rho),5)))
/(this->Tcrit()*pow((pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)),3)*pow(((this->W()/rho)-this->b()),3));
}




//(molar values)
//using second order implicit differentiation 
inline scalar pengRobinson::d2vdT2(const scalar rho, const scalar T) const
{
return -(pow(this->dpdT(rho,T),2)*this->d2pdv2(rho,T)+   pow(this->dpdv(rho,T),2)  *this->d2pdT2(rho,T)- 2*this->dpdv(rho,T)*this->dpdT(rho,T)*this->d2pdvdT(rho,T))/(   pow(this->dpdv(rho,T),3));
}




//(molar values)
inline scalar pengRobinson::d2pdvdT(const scalar rho, const scalar T) const
{
return -2*this->a()*this->n()*(this->b()+(this->W()/rho))*(this->n()+1)*pow(T/this->Tcrit(),0.5)
/(T*pow((pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2)),2))+(2*this->a()*pow(this->n(),2)*(this->b()+(this->W()/rho)))/(this->Tcrit()*pow(pow(this->b(),2)-2*this->b()*(this->W()/rho)-pow((this->W()/rho),2),2))-this->RR/(pow(this->b()-(this->W()/rho),2));
}




// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar pengRobinson::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{
return pow(2,0.5)*this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)
*log(this->b()*(pow(2,0.5)+1)+(this->W()/rho))/(8*this->b()*pow(T,2))
-pow(2,0.5)*this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)
*log(this->b()*(1-pow(2,0.5))+(this->W()/rho))/(8*this->b()*pow(T,2));
}





//* * * * * * * * * * * * * thermodynamic properties * * * * * * * * * * * //


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar pengRobinson::isobarExpCoef(const scalar rho,const scalar T)  const
{
 return this->dvdT(rho, T)*rho/this->W();
}

//isothemal compressiblity kappa 
//(molar values)
inline scalar pengRobinson::isothermalCompressiblity(const scalar rho,const scalar T) const
{
return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
//also possible : return -this->dvdp(rho,T)*rho/this->W();
}










// Construct from components
inline pengRobinson::pengRobinson
(
    const specie& sp,
    scalar pcrit, 
    scalar Tcrit,
    scalar azentricFactor
    
)
:
    specie(sp),    
    pcrit_(pcrit),
    Tcrit_(Tcrit),
    azentricFactor_(azentricFactor)
{
	rhostd_=this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R()));	// Starting GUESS for the density by ideal gas law
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy

inline pengRobinson::pengRobinson(const word& name,  pengRobinson& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
    azentricFactor_(pg.azentricFactor_)
{
  pg.rhostd_=this->rho(Pstd,Tstd,  (Pstd*this->W()/(Tstd*this->R())));	// Starting GUESS for the density by ideal gas law
}


// Construct and return a clone
inline autoPtr<pengRobinson> pengRobinson::clone() const
{
    return autoPtr<pengRobinson>(new pengRobinson(*this));
}


// Selector from Istream
inline autoPtr<pengRobinson> pengRobinson::New(Istream& is)
{
    return autoPtr<pengRobinson>(new pengRobinson(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return density [kg/m^3]
// TO DO Include a max Iteration number loop and abort function 
inline scalar pengRobinson::rho(const scalar p,const scalar T,const scalar rho0) const
{

   scalar molarVolumePrevIteration;
   scalar molarVolume;
   int iter=0;	
   int i;
   int maxIter_=400;
   scalar tol_=1e-8;
   scalar rho1=rhoMax_, rho2=rhoMin_,rho3, f1,f2,f3;

molarVolume=this->W()/rho0;	 

   do
   {
   molarVolumePrevIteration= molarVolume;
   i=0;
      	do
   	{
        
  	molarVolume=molarVolumePrevIteration-((this->pReturn((this->W()/molarVolumePrevIteration),T)-p)/(this->dpdv((this->W()/       	
  	molarVolumePrevIteration),T)))/(pow(2,i)); 


	i++;
		if(i>8) 
		{
		  //using bisection methode as backup, solution must be between rho=0.001 to rho=1500; 
		  for(i=0;i<200;i++)
		  {
		
                  	f1= (this->pReturn(rho1,T)-p);
                  	f2= (this->pReturn(rho2,T)-p);
                  	rho3=(rho1+rho2)/2;
		  	f3=(this->pReturn(rho3,T)-p);
	
			if ((f2<0 && f3>0)||(f2>0 &&f3<0))	
		  	{
		 		rho1=rho3;
			}  
 		  	else if ((f1<0 && f3>0)||(f1>0 &&f3<0))
                  	{
		 		rho2=rho3;
			}
		 
			else	
			{
                        	rho2=(rho2+rho3)/2;
			}

			if(mag(f3)<p*tol_)
                  	{	  	  
				molarVolume=this->W()/rho3;
                 	 	molarVolumePrevIteration=this->W()/rho3;
				break; 
		  	}
			else
			{
				molarVolumePrevIteration=this->W()/rho3;
			}
		    }		  
  		}


        }
        while(mag(this->pReturn((this->W()/molarVolume),T)-p)>mag(this->pReturn((this->W()/molarVolumePrevIteration),T)-p));

           if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "inline scalar pengRobinson::rho(const scalar p,const scalar T,const scalar rho0) const "
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

   }
   while(mag(molarVolumePrevIteration-molarVolume)>tol_*(this->W()/rho0));





   return this->W()/molarVolume;   
 
}

//- Return density [kg/m^3]on 
inline scalar pengRobinson::rho(const scalar p,const scalar T) const
{
	scalar rho0=p/(this->R()*T); //using perfect gas equation as starting point
	return rho(p,T,rho0); 
}

//- Return compressibility drho/dp [s^2/m^2] 
inline scalar pengRobinson::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho,T)*pow(rho,2)/this->W();
}

//- Return compression factor []
inline scalar pengRobinson::Z( const scalar p, const scalar T,const scalar rho0) const
{
    return (p*this->rho(p,T,rho0))/(this->R()*T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void pengRobinson::operator+=(const pengRobinson& pg)
{
    specie::operator+=(pg);
}

	inline void pengRobinson::operator-=(const pengRobinson& pg)
{
    specie::operator-=(pg);
}

inline void pengRobinson::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
//****************not working**************//
/*
inline pengRobinson operator+
(
    const pengRobinson& pg1,
    const pengRobinson& pg2
)
{
    return pengRobinson
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline pengRobinson operator-
(
    const pengRobinson& pg1,
    const pengRobinson& pg2
)
{
    return pengRobinson
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline pengRobinson operator*
(
    const scalar s,
    const pengRobinson& pg
)
{
    return pengRobinson(s*static_cast<const specie&>(pg));
}


inline pengRobinson operator==
(
    const pengRobinson& pg1,
    const pengRobinson& pg2
)
{
    return pg2 - pg1;
}

*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
