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

#include "aungierRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline aungierRedlichKwong::aungierRedlichKwong
(
    const specie& sp,
    scalar pcrit, 
    scalar Tcrit,
    scalar azentricFactor,
    scalar rhocrit
)
:
    specie(sp),    
    pcrit_(pcrit),
    Tcrit_(Tcrit),
    azentricFactor_(azentricFactor),
    rhocrit_(rhocrit),
    a_(0.42747*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
    c_(this->RR*Tcrit_/(pcrit_+(a_/(this->W()/rhocrit_*(this->W()/rhocrit_+b_))))+b_-this->W()/rhocrit_),
    n_(0.4986+1.2735*azentricFactor_+0.4754*pow(azentricFactor_,2)),
    // Starting GUESS for the density by ideal gas law
        rhostd_(this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R())))
{}


// Construct as named copy
inline aungierRedlichKwong::aungierRedlichKwong(const word& name, const aungierRedlichKwong& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
    azentricFactor_(pg.azentricFactor_),
    rhocrit_(pg.rhocrit_),    
    a_(pg.a_),
    b_(pg.b_),
    c_(pg.c_),
    n_(pg.n_),
    rhostd_(pg.rhostd_)
{}



// Construct and return a clone
inline autoPtr<aungierRedlichKwong> aungierRedlichKwong::clone() const
{
    return autoPtr<aungierRedlichKwong>(new aungierRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<aungierRedlichKwong> aungierRedlichKwong::New(Istream& is)
{
    return autoPtr<aungierRedlichKwong>(new aungierRedlichKwong(is));
}


// * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * //
inline scalar aungierRedlichKwong::rhostd()const
{
    return rhostd_;
}

//returns the pressure for a given density and temperature    
inline scalar aungierRedlichKwong::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return this->RR*T/(Vm-b_+c_)
        -a_*pow(T,-n_)/(pow(Tcrit_,-n_)*Vm*(Vm+b_));
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar aungierRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return  a_*pow(T,-n_)*pow(Tcrit_,n_)*(b_+2*Vm)/(pow(Vm,2)*pow((b_+Vm),2))
        -this->RR*T/(pow((b_-Vm-c_),2));
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar aungierRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return a_*n_*pow(T,-n_-1)*pow(Tcrit_,n_)/(Vm*(Vm+b_))
        -this->RR/(b_-Vm-c_);
}


//Real  deviative dv/dT at constant pressure 
//using implicit differentiation
//(molar values)
inline scalar aungierRedlichKwong::dvdT(const scalar rho,const scalar T) const
{
    return (-1)*this->dpdT(rho,T)/this->dpdv(rho,T);
}


//Real  deviative dv/dp at constant temperature 
//(molar values)
inline scalar aungierRedlichKwong::dvdp(const scalar rho,const scalar T)  const 
{
    return 1/this->dpdv(rho,T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar aungierRedlichKwong::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    return -pow((T/Tcrit_),-n_)*(a_*log(Vm)/(b_)
        -a_*log(b_+Vm)/b_)+this->RR*T*log(Vm-b_+c_);
}



//needed to calculate the entropy
//(molar values)
inline scalar aungierRedlichKwong::integral_dpdT_dv
( 
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    return pow(T,-n_-1)*pow(Tcrit_,n_)*(a_*n_*log(Vm)/(b_) 
        -a_*n_*log(b_+Vm)/(b_))+ this->RR*log(-b_+c_+Vm);
}


//(molar values)
inline scalar aungierRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -a_*n_*pow(T,-n_-2)*pow(Tcrit_,n_)*(n_+1)/(Vm*(b_+Vm));
}




//(molar values)
inline scalar aungierRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -2*a_*pow(T,-n_)*pow(Tcrit_,n_)*(pow(b_,2)
        +3*b_*Vm+3*pow(Vm,2))/(pow(Vm,3)*pow((b_+Vm),3))
        -2*this->RR*T/(pow((b_-Vm-c_),3));
}




//(molar values)
//using second order implicit differentiation 
inline scalar aungierRedlichKwong::d2vdT2(const scalar rho, const scalar T) const
{
    return 
    -(
          pow(this->dpdT(rho,T),2)*this->d2pdv2(rho,T)
        + pow(this->dpdv(rho,T),2)*this->d2pdT2(rho,T)
        - 2*this->dpdv(rho,T)*this->dpdT(rho,T)*this->d2pdvdT(rho,T)
    )
    /(pow(this->dpdv(rho,T),3));
}




//(molar values)
inline scalar aungierRedlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -a_*n_*pow(T,-n_-1)*pow(Tcrit_,n_)*(b_+2*Vm)/(pow(Vm,2)*pow((b_+Vm),2))
        -this->RR/(pow((b_-c_-Vm),2));
}




// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar aungierRedlichKwong::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
) const 
{
    scalar Vm = this->W()/rho;
    return pow(T,-n_-2)*pow(Tcrit_,n_)*(a_*n_*(n_+1)*log(b_+Vm)
        /b_-a_*n_*(1+n_)*log(Vm)/b_);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar aungierRedlichKwong::isobarExpCoef(const scalar rho,const scalar T)  const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar aungierRedlichKwong::isothermalCompressiblity(const scalar rho,const scalar T) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //also possible : return -this->dvdp(rho,T)*rho/this->W();
}


//- Return density [kg/m^3]
inline scalar aungierRedlichKwong::rho(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
 
   scalar molarVolumePrevIteration;
   scalar molarVolume;
   int iter=0;
   int maxIter_=400;
   scalar tol_=1e-8;
   scalar rho1=rhoMax_;
   scalar rho2=rhoMin_;

   molarVolume=this->W()/rho0;	 

   do
   {
       molarVolumePrevIteration= molarVolume;
   
       label i=0;
       do
       {
           molarVolume=molarVolumePrevIteration
               -(
                   (this->p((this->W()/molarVolumePrevIteration),T) - p)
                   /(this->dpdv((this->W()/molarVolumePrevIteration),T))
               )/pow(2,i);

           i++;
           if (i>8) 
           {
               //using bisection methode as backup,
               //solution must be between rho=0.001 to rho=1500;
               for(i=0; i<200; i++)
               {
                   scalar f1 = this->p(rho1,T) - p;
                   scalar f2 = this->p(rho2,T) - p;
                   scalar rho3 = (rho1 + rho2)/2;
                   scalar f3 = this->p(rho3,T) - p;
	
                   if ((f2 < 0 && f3 > 0) || (f2 > 0 && f3 < 0))
                   {
                       rho1=rho3;
                   }  
                   else if ((f1 < 0 && f3 > 0)||(f1 > 0 && f3 < 0))
                   {
                       rho2=rho3;
                   }		 
                   else	
                   {
                       rho2=(rho2 + rho3)/2;
                   }

                   if(mag(f3) < p*tol_)
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
       while
       (
           mag(this->p((this->W()/molarVolume),T) - p)
         > mag(this->p((this->W()/molarVolumePrevIteration),T) - p)
       );

       if (iter++ > maxIter_)
       {
           FatalErrorIn
           (
               "inline scalar redlichKwong::rho(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

   return this->W()/molarVolume;
}

//- Return density [kg/m^3]
inline scalar aungierRedlichKwong::rho(const scalar p,const scalar T) const
{
    //using perfect gas equation as starting point
    return rho(p,T,p/(this->R()*T));
}

//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar aungierRedlichKwong::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho,T)*pow(rho,2)/this->W();
}

//- Return compression factor []
inline scalar aungierRedlichKwong::Z( const scalar p, const scalar T,const scalar rho0) const
{
    return (p*this->rho(p,T,rho0))/(this->R()*T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void aungierRedlichKwong::operator+=(const aungierRedlichKwong& pg)
{
    specie::operator+=(pg);
}

	inline void aungierRedlichKwong::operator-=(const aungierRedlichKwong& pg)
{
    specie::operator-=(pg);
}

inline void aungierRedlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
//****************not working**************//
/*
inline aungierRedlichKwong operator+
(
    const aungierRedlichKwong& pg1,
    const aungierRedlichKwong& pg2
)
{
    return aungierRedlichKwong
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline aungierRedlichKwong operator-
(
    const aungierRedlichKwong& pg1,
    const aungierRedlichKwong& pg2
)
{
    return aungierRedlichKwong
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline aungierRedlichKwong operator*
(
    const scalar s,
    const aungierRedlichKwong& pg
)
{
    return aungierRedlichKwong(s*static_cast<const specie&>(pg));
}


inline aungierRedlichKwong operator==
(
    const aungierRedlichKwong& pg1,
    const aungierRedlichKwong& pg2
)
{
    return pg2 - pg1;
}

*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
