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

#include "soaveRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{





// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline scalar soaveRedlichKwong::pcrit()const
{
return pcrit_;
}



inline scalar soaveRedlichKwong::Tcrit()const
{
return Tcrit_;
}


inline scalar soaveRedlichKwong::rhostd()const
{
return rhostd_;
}


// Returns the Azentric Factor (Acentric Factor)
inline scalar soaveRedlichKwong::azentricFactor() const
{
return azentricFactor_;
}

//returns the pressure for a given density and temperature    
inline scalar soaveRedlichKwong::pReturn(const scalar rho,const scalar T) const
{
return  (this->RR*T/((this->W()/rho)-this->b())-this->a()*pow((1+this->n()*(1-pow((T/this->Tcrit()),0.5))),2)/((this->W()/rho)*((this->W()/rho)+this->b())));
}



// Faktor a of the soave redlich Kwong equation of state  
//(molar values)
inline scalar soaveRedlichKwong::a() const
{
return 0.42747*pow(this->RR,2)*pow(this->Tcrit(),2)/(this->pcrit());
}



// Faktor b of the soave redlich Kwong equation of state 
//(molar values)
inline scalar soaveRedlichKwong::b() const
{
return 0.08664*this->RR*this->Tcrit()/this->pcrit(); 
}

// Faktor n of the soave redlich Kwong equation of state  
//(molar values)
inline scalar soaveRedlichKwong::n() const
{
return 0.48508+1.55171*this->azentricFactor()-0.15613*pow(this->azentricFactor(),2);
}


//* * * * * * * * * * * * * Derivatives * * * * * * * * * * * //



//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar soaveRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
 return -(2*this->a()*this->n()*this->Tcrit()*(this->b()-(this->W()/rho))*(pow(this->b(),2)+this->b()*(this->W()/rho)-2*pow((this->W()/rho),2))*(this->n()+1)*pow((T/this->Tcrit()),0.5)+this->Tcrit()*(this->RR*T*pow((this->W()/rho),2)*(pow(this->b(),2)+2*this->b()*(this->W()/rho)+pow((this->W()/rho),2))-this->a()*(pow(this->b(),3)-3*this->b()*pow((this->W()/rho),2)+2*pow((this->W()/rho),3))*pow((this->n()+1),2))-this->a()*pow(this->n(),2)*T*(pow(this->b(),3)-3*this->b()*pow((this->W()/rho),2)+2*pow((this->W()/rho),3)))/(pow((this->W()/rho),2)*this->Tcrit()*pow((this->b()+(this->W()/rho)),2)*pow((this->b()-(this->W()/rho)),2));
}




//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar soaveRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
return this->a()*this->n()*(this->n()+1)*pow((T/this->Tcrit()),0.5)/(T*(this->W()/rho)*(this->b()+(this->W()/rho)))-this->a()*pow(this->n(),2)/((this->W()/rho)*this->Tcrit()*(this->b()+(this->W()/rho)))-this->RR/(this->b()-(this->W()/rho));
}




//Real  deviative dv/dT at constant pressure
//using implicit differentiation
// (molar values)
inline scalar soaveRedlichKwong::dvdT(const scalar rho,const scalar T) const
{
return (-1)*this->dpdT(rho,T)/this->dpdv(rho,T);
}




//Real  deviative dv/dp at constant temperature 
//(molar values)
inline scalar soaveRedlichKwong::dvdp(const scalar rho,const scalar T)  const 
{
return 1/this->dpdv(rho,T);
}



//needed to calculate the internal energy 
//(molar values)
inline scalar soaveRedlichKwong::integral_p_dv(const scalar rho,const scalar T) const
{
return this->RR*T*log((this->W()/rho)-this->b())
-(this->a()*(2*this->n()*this->Tcrit()*(this->n()+1)*pow(T/this->Tcrit(),0.5)-this->Tcrit()*(pow(this->n(),2)+2*this->n()+1)-pow(this->n(),2)*T)*log(this->b()+(this->W()/rho)))/(this->b()*this->Tcrit())
+this->a()*(2*this->n()*this->Tcrit()*(this->n()+1)*pow(T/this->Tcrit(),0.5)-this->Tcrit()*(pow(this->n(),2)+2*this->n()+1)-pow(this->n(),2)*T)*log((this->W()/rho))/(this->b()*this->Tcrit());
}



//needed to calculate the entropy
//(molar values)
inline scalar soaveRedlichKwong::integral_dpdT_dv(const scalar rho,const scalar T) const
{
return this->RR*log((this->W()/rho)-this->b())+(pow(this->n(),2)*this->a()/(this->b()*this->Tcrit())-this->a()*this->n()*(this->n()+1)*pow((T/this->Tcrit()),0.5)/(this->b()*T))*log(this->b()+(this->W()/rho))+(this->a()*this->n()*(this->n()+1)*pow((T/this->Tcrit()),0.5)/(this->b()*T)-this->a()*pow(this->n(),2)/(this->b()*this->Tcrit()))*log((this->W()/rho));
}




//* * * * * * * * * * * * * second order Derivative based functions * * * * * * * * * * * //



//(molar values)
inline scalar soaveRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
return -this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)/(2*pow(T,2)*(this->W()/rho)*(this->b()+(this->W()/rho)));
}




//(molar values) 
inline scalar soaveRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
 return 2*(2*this->a()*this->n()*this->Tcrit()*(this->b()-(this->W()/rho))
    *(pow(this->b(),4)+pow(this->b(),3)*(this->W()/rho)-2*pow(this->b(),2)*pow((this->W()/rho),2)
	-3*this->b()*pow((this->W()/rho),3)+3*pow((this->W()/rho),4))*(this->n()+1)*pow(T/this->Tcrit(),0.5)
	-this->Tcrit()*(this->a()*(pow(this->b(),5)-3*pow(this->b(),3)*pow((this->W()/rho),2)-pow(this->b(),2)*pow((this->W()/rho),3)
	+6*this->b()*pow((this->W()/rho),4)-3*pow((this->W()/rho),5))*pow((this->n()+1),2)
	+this->RR*T*pow((this->W()/rho),3)*(pow(this->b(),3)+3*pow(this->b(),2)*(this->W()/rho)+3*this->b()*pow((this->W()/rho),2)+pow((this->W()/rho),3)))
	-this->a()*pow(this->n(),2)*T*(pow(this->b(),5)-3*pow(this->b(),3)*pow((this->W()/rho),2)-pow(this->b(),2)*pow((this->W()/rho),3)
	+6*this->b()*pow((this->W()/rho),4)-3*pow((this->W()/rho),5)))
	/(pow((this->W()/rho),3)*this->Tcrit()*pow(((this->W()/rho)+this->b()),3)*pow((this->b()-(this->W()/rho)),3));
}




//(molar values)
// using second Order implicit differentiation 
inline scalar soaveRedlichKwong::d2vdT2(const scalar rho, const scalar T) const
{
return -(pow(this->dpdT(rho,T),2)*this->d2pdv2(rho,T)+   pow(this->dpdv(rho,T),2)  *this->d2pdT2(rho,T)- 2*this->dpdv(rho,T)*this->dpdT(rho,T)*this->d2pdvdT(rho,T))/(   pow(this->dpdv(rho,T),3));
}




//(molar values)
inline scalar soaveRedlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
return -this->a()*this->n()*(this->b()+2*(this->W()/rho))*(this->n()+1)*pow(T/this->Tcrit(),0.5)
/(T*pow((this->W()/rho),2)*pow((this->b()+(this->W()/rho)),2))
-(this->RR*pow((this->W()/rho),2)*this->Tcrit()*(pow(this->b(),2)+2*this->b()*(this->W()/rho)+pow((this->W()/rho),2))
  -this->a()*pow(this->n(),2)*(pow(this->b(),3)-3*this->b()*pow((this->W()/rho),2)+2*pow((this->W()/rho),3)))
  /(pow((this->W()/rho),2)*this->Tcrit()*pow((this->b()+(this->W()/rho)),2)*pow((this->b()-(this->W()/rho)),2));
}




// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar soaveRedlichKwong::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{
return this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)*log(this->b()+(this->W()/rho))
/(2*this->b()*pow(T,2))
-this->a()*this->n()*(this->n()+1)*pow(T/this->Tcrit(),0.5)*log((this->W()/rho))/(2*this->b()*pow(T,2));
}





//* * * * * * * * * * * * * thermodynamic properties * * * * * * * * * * * //


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar soaveRedlichKwong::isobarExpCoef(const scalar rho,const scalar T)  const
{
return this->dvdT(rho, T)*rho/this->W();
}

//isothemal compressiblity kappa 
//(molar values)
inline scalar soaveRedlichKwong::isothermalCompressiblity(const scalar rho,const scalar T) const
{
return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
}










// Construct from components
inline soaveRedlichKwong::soaveRedlichKwong
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
	rhostd_=this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy
inline soaveRedlichKwong::soaveRedlichKwong(const word& name,  soaveRedlichKwong& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
	azentricFactor_(pg.azentricFactor_)
{
  pg.rhostd_=this->rho(Pstd,Tstd,  (Pstd*this->W()/(Tstd*this->R())));	
}


// Construct and return a clone
inline autoPtr<soaveRedlichKwong> soaveRedlichKwong::clone() const
{
    return autoPtr<soaveRedlichKwong>(new soaveRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<soaveRedlichKwong> soaveRedlichKwong::New(Istream& is)
{
    return autoPtr<soaveRedlichKwong>(new soaveRedlichKwong(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return density [kg/m^3]on 
inline scalar soaveRedlichKwong::rho(const scalar p,const scalar T,const scalar rho0) const
{
 
   scalar molarVolumePrevIteration;
   scalar molarVolume;
   int iter=0;	
   int maxIter_=400;
   scalar tol_=1e-8;
   int i;
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
               "inline scalar soaveRedlichKwong::rho(const scalar p,const scalar T,const scalar rho0) const "
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

   }
   while(mag(molarVolumePrevIteration-molarVolume)>tol_*(this->W()/rho0));

   return this->W()/molarVolume;   

}

//- Return density [kg/m^3]on 
inline scalar soaveRedlichKwong::rho(const scalar p,const scalar T) const
{

	scalar rho0=p/(this->R()*T); //using perfect gas equation as starting point
	return rho(p,T,rho0);
}

//- Return compressibility drho/dp [s^2/m^2]
inline scalar soaveRedlichKwong::psi(const scalar rho, const scalar T) const
{
  return -this->dvdp(rho,T)*pow(rho,2)/this->W();
}

//- Return compression factor []
inline scalar soaveRedlichKwong::Z( const scalar p, const scalar T,const scalar rho0) const
{
        return (p*this->rho(p,T,rho0))/(this->R()*T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void soaveRedlichKwong::operator+=(const soaveRedlichKwong& pg)
{
    specie::operator+=(pg);
}

	inline void soaveRedlichKwong::operator-=(const soaveRedlichKwong& pg)
{
    specie::operator-=(pg);
}

inline void soaveRedlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
//**************** not working**************//
/*
inline soaveRedlichKwong operator+
(
    const soaveRedlichKwong& pg1,
    const soaveRedlichKwong& pg2
)
{
    return soaveRedlichKwong
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline soaveRedlichKwong operator-
(
    const soaveRedlichKwong& pg1,
    const soaveRedlichKwong& pg2
)
{
    return soaveRedlichKwong
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline soaveRedlichKwong operator*
(
    const scalar s,
    const soaveRedlichKwong& pg
)
{
    return soaveRedlichKwong(s*static_cast<const specie&>(pg));
}


inline soaveRedlichKwong operator==
(
    const soaveRedlichKwong& pg1,
    const soaveRedlichKwong& pg2
)
{
    return pg2 - pg1;
}

*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
