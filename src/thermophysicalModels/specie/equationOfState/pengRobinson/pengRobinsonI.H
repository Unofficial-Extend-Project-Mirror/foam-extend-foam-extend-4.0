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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    azentricFactor_(azentricFactor),
    a_(0.457235*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_),
    b_(0.077796*this->RR*Tcrit_/pcrit_), 
    n_(0.37464+1.54226*azentricFactor_-0.26992*pow(azentricFactor_,2)), 
    rhostd_(this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R()))) 
{}

// Construct as named copy
inline pengRobinson::pengRobinson(const word& name, const pengRobinson& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
    azentricFactor_(pg.azentricFactor_),
    a_(pg.a_),
    b_(pg.b_),
    n_(pg.n_),
    rhostd_(pg.rhostd_)
{}


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


// * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * //

inline scalar pengRobinson::rhostd()const
{
    return rhostd_;
}


//returns the pressure for a given density and temperature    
inline scalar pengRobinson::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return  this->RR*T/(Vm-b_)
        -(a_*pow((1+n_*(1-pow((T/Tcrit_),0.5))),2))
        /(pow(Vm,2)+2*b_*Vm-pow(b_,2));
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar pengRobinson::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    -(
        4*a_*n_*Tcrit_*(b_-Vm)*(pow(b_,2)
        -pow(Vm,2))*(n_+1)*pow((T/Tcrit_),0.5)
        +Tcrit_*(-2*a_*pow((n_+1),2)*(pow(b_,3)-pow(b_,2)*Vm
        -b_*pow(Vm,2)+pow(Vm,3))
        +this->RR*T*(pow(b_,4)-4*pow(b_,3)*Vm
        +2*pow(b_,2)*pow(Vm,2)
        +4*b_*pow(Vm,3)+pow(Vm,4)))
        -2*a_*pow(n_,2)*T*(pow(b_,3)-pow(b_,2)*Vm
        -b_*pow(Vm,2)+pow(Vm,3))
     )
     /(Tcrit_*pow((b_-Vm),2)*pow((pow(b_,2)
        -2*b_*Vm-pow(Vm,2)),2));
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar pengRobinson::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    (
        -a_*n_*(n_+1)*pow((T/Tcrit_),0.5)
            /(T*(pow(b_,2)
        -2*b_*Vm-pow(Vm,2)))
        +a_*pow(n_,2)/(Tcrit_*(pow(b_,2)
        -2*b_*Vm-pow(Vm,2)))
        -this->RR/(b_-Vm)
    );
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
    scalar Vm = this->W()/rho;
    return  pow(2,0.5)*a_*(2*n_*Tcrit_*(n_+1)*pow(T/Tcrit_,0.5)
        -Tcrit_*(pow(n_,2)+2*n_+1)-pow(n_,2)*T)
            *log(b_*(1-pow(2,0.5))+Vm)/(4*b_*Tcrit_)
        +this->RR*T*log(Vm-b_)
        -pow(2,0.5)*a_*(2*n_*Tcrit_*(n_+1)*pow(T/Tcrit_,0.5)
        -Tcrit_*(pow(n_,2)+2*n_+1)
        -pow(n_,2)*T)*log(b_*(pow(2,0.5)+1)+Vm)
            /(4*b_*Tcrit_);
}


//(molar values)
//needed to calculate the entropy
inline scalar pengRobinson::integral_dpdT_dv(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return (pow(2,0.5)*a_*n_*(n_+1)*pow(T/Tcrit_,0.5)/(4*b_*T)
        -pow(2,0.5)*a_*pow(n_,2)/(4*b_*Tcrit_))
            *log(b_*(1-pow(2,0.5))+Vm)
        +this->RR*log(Vm-b_)
        +(pow(2,0.5)*a_*pow(n_,2)/(4*b_*Tcrit_)
        -pow(2,0.5)*a_*n_*(n_+1)*pow(T/Tcrit_,0.5)/(4*b_*T))
            *log(b_*(pow(2,0.5)+1)+Vm);
}


//(molar values)
inline scalar pengRobinson::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    a_*n_*(n_+1)*pow((T/Tcrit_),0.5)
        /(2*pow(T,2)*(pow(b_,2)
    -2*b_*Vm-pow(Vm,2)));
}


//(molar values)
inline scalar pengRobinson::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    2*(
        2*a_*n_*Tcrit_*(b_-Vm)*(n_+1)*
        (
            5*pow(b_,4)-4*pow(b_,3)*Vm
            -4*pow(b_,2)*pow(Vm,2)
            +3*pow(Vm,4)
        )
        *pow(T/Tcrit_,0.5)
        +Tcrit_*
        (
            this->RR*T*
            (
                pow(b_,6)
                -6*pow(b_,5)*Vm
                +9*pow(b_,4)*pow(Vm,2)
                +4*pow(b_,3)*pow(Vm,3)
                -9*pow(b_,2)*pow(Vm,4)
                -6*b_*pow(Vm,5)
                -pow(Vm,6)
            )
            -a_*pow((n_+1),2)*
            (
                 5*pow(b_,5)
                 -9*pow(b_,4)*Vm
                 +4*pow(b_,2)*pow(Vm,3)
                 +3*b_*pow(Vm,4)-3*pow(Vm,5)
            )
         )
         -a_*pow(n_,2)*T*
         (
             5*pow(b_,5)
             -9*pow(b_,4)*Vm
             +4*pow(b_,2)*pow(Vm,3)
             +3*b_*pow(Vm,4)
             -3*pow(Vm,5)
         )
     )
     /
     (
         Tcrit_*
         pow((pow(b_,2)
             -2*b_*Vm-pow(Vm,2)),3)
         *pow((Vm-b_),3)
     );
}


//(molar values)
//using second order implicit differentiation 
inline scalar pengRobinson::d2vdT2(const scalar rho, const scalar T) const
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
inline scalar pengRobinson::d2pdvdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -2*a_*n_*(b_+Vm)*(n_+1)*pow(T/Tcrit_,0.5)
        /(T*pow((pow(b_,2)-2*b_*Vm-pow(Vm,2)),2))
    +(2*a_*pow(n_,2)*(b_+Vm))
        /(Tcrit_*pow(pow(b_,2)-2*b_*Vm-pow(Vm,2),2))
    -this->RR/(pow(b_-Vm,2));
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar pengRobinson::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{
    scalar Vm = this->W()/rho;
    return pow(2,0.5)*a_*n_*(n_+1)*pow(T/Tcrit_,0.5)
        *log(b_*(pow(2,0.5)+1)+Vm)/(8*b_*pow(T,2))
    -pow(2,0.5)*a_*n_*(n_+1)*pow(T/Tcrit_,0.5)
        *log(b_*(1-pow(2,0.5))+Vm)/(8*b_*pow(T,2));
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar pengRobinson::isobarExpCoef(const scalar rho,const scalar T)  const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar pengRobinson::isothermalCompressiblity(const scalar rho,const scalar T) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //also possible : return -this->dvdp(rho,T)*rho/this->W();
}


//- Return density [kg/m^3]
inline scalar pengRobinson::rho(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
 
   scalar molarVolumePrevIteration;
   scalar molarVolume;
   label iter=0;
   label maxIter_=400;
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


//- Return density [kg/m^3]on 
inline scalar pengRobinson::rho(const scalar p,const scalar T) const
{
    // using perfect gas equation as starting point
    return rho(p,T,p/(this->R()*T)); 
}


//- Return compressibility drho/dp at T=constant [s^2/m^2] 
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
