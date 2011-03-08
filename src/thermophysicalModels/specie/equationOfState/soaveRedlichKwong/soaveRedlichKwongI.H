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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    azentricFactor_(azentricFactor),
    a_(0.42747*pow(this->RR,2)*pow(Tcrit_,2)/(pcrit_)),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
    n_(0.48508+1.55171*azentricFactor_-0.15613*pow(azentricFactor_,2)),
        // Starting GUESS for the density by ideal gas law
        rhostd_(this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R())))
{}




// Construct as named copy
inline soaveRedlichKwong::soaveRedlichKwong(const word& name,const  soaveRedlichKwong& pg)
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
inline autoPtr<soaveRedlichKwong> soaveRedlichKwong::clone() const
{
    return autoPtr<soaveRedlichKwong>(new soaveRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<soaveRedlichKwong> soaveRedlichKwong::New(Istream& is)
{
    return autoPtr<soaveRedlichKwong>(new soaveRedlichKwong(is));
}


// * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * //

inline scalar soaveRedlichKwong::rhostd()const
{
    return rhostd_;
}


//returns the pressure for a given density and temperature    
inline scalar soaveRedlichKwong::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    (
        this->RR*T/(Vm-b_)
        -a_*pow((1+n_*(1-pow((T/Tcrit_),0.5))),2)
            /(Vm*(Vm+b_))
    );
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar soaveRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    -(
        2*a_*n_*Tcrit_*(b_-Vm)*(pow(b_,2)
        +b_*Vm-2*pow(Vm,2))*(n_+1)*pow((T/Tcrit_),0.5)
        +Tcrit_*
        (
            this->RR*T*pow(Vm,2)*(pow(b_,2)
                +2*b_*Vm+pow(Vm,2))
        -a_*
            (pow(b_,3)-3*b_*pow(Vm,2)+
            2*pow(Vm,3))
        *pow((n_+1),2)
        )
        -a_*pow(n_,2)*T*(pow(b_,3)
        -3*b_*pow(Vm,2)+2*pow(Vm,3))
     )
     /
     (
         pow(Vm,2)*Tcrit_*pow((b_+Vm),2)
             *pow((b_-Vm),2)
     );
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar soaveRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return a_*n_*(n_+1)*pow((T/Tcrit_),0.5)
        /(T*Vm*(b_+Vm))
    -a_*pow(n_,2)
        /(Vm*Tcrit_*(b_+Vm))
    -this->RR/(b_-Vm);
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
    scalar Vm = this->W()/rho;
    return this->RR*T*log(Vm-b_)
    -(
        a_*(2*n_*Tcrit_*(n_+1)*pow(T/Tcrit_,0.5)
        -Tcrit_*(pow(n_,2)+2*n_+1)
        -pow(n_,2)*T)*log(b_+Vm)
     )
     /(b_*Tcrit_)
     +a_*
     (
         2*n_*Tcrit_*(n_+1)*pow(T/Tcrit_,0.5)
         -Tcrit_*(pow(n_,2)+2*n_+1)-pow(n_,2)*T
     )
     *log(Vm)/(b_*Tcrit_);
}


//needed to calculate the entropy
//(molar values)
inline scalar soaveRedlichKwong::integral_dpdT_dv(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return this->RR*log(Vm-b_)
    +(
        pow(n_,2)*a_/(b_*Tcrit_)
        -a_*n_*(n_+1)*pow((T/Tcrit_),0.5)
            /(b_*T)
    )*log(b_+Vm)
    +(a_*n_*(n_+1)*pow((T/Tcrit_),0.5)
        /(b_*T)
    -a_*pow(n_,2)/(b_*Tcrit_))*log(Vm);
}


//(molar values)
inline scalar soaveRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -a_*n_*(n_+1)*pow(T/Tcrit_,0.5)
        /(2*pow(T,2)*Vm*(b_+Vm));
}


//(molar values) 
inline scalar soaveRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 
    2*(
        2*a_*n_*Tcrit_*(b_-Vm)*
        (
            pow(b_,4)+pow(b_,3)*Vm
            -2*pow(b_,2)*pow(Vm,2)
	    -3*b_*pow(Vm,3)
            +3*pow(Vm,4)
        )
        *(n_+1)*pow(T/Tcrit_,0.5)
	-Tcrit_*
        (
            a_*
            (
                pow(b_,5)
                -3*pow(b_,3)*pow(Vm,2)
                -pow(b_,2)*pow(Vm,3)
	        +6*b_*pow(Vm,4)
                -3*pow(Vm,5)
            )*pow((n_+1),2)
	    +this->RR*T*pow(Vm,3)*
            (
                pow(b_,3)
                +3*pow(b_,2)*Vm
                +3*b_*pow(Vm,2)
                +pow(Vm,3)
            )
        )
	-a_*pow(n_,2)*T*
        (
            pow(b_,5)
            -3*pow(b_,3)*pow(Vm,2)
            -pow(b_,2)*pow(Vm,3)
	    +6*b_*pow(Vm,4)
            -3*pow(Vm,5)
        )
    )
    /(pow(Vm,3)*Tcrit_*pow((Vm+b_),3)
        *pow((b_-Vm),3));
}


//(molar values)
// using second Order implicit differentiation 
inline scalar soaveRedlichKwong::d2vdT2(const scalar rho, const scalar T) const
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
inline scalar soaveRedlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -a_*n_*(b_+2*Vm)*(n_+1)*pow(T/Tcrit_,0.5)
        /(T*pow(Vm,2)*pow((b_+Vm),2))
    -(
        this->RR*pow(Vm,2)*Tcrit_*(pow(b_,2)
        +2*b_*Vm+pow(Vm,2))
        -a_*pow(n_,2)*
        (
            pow(b_,3)
            -3*b_*pow(Vm,2)
            +2*pow(Vm,3)
        )
    )
  /(pow(Vm,2)*Tcrit_*pow((b_+Vm),2)*pow((b_-Vm),2));
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar soaveRedlichKwong::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{    
    scalar Vm = this->W()/rho;
    return a_*n_*(n_+1)*pow(T/Tcrit_,0.5)*log(b_+Vm)
        /(2*b_*pow(T,2))
    -a_*n_*(n_+1)*pow(T/Tcrit_,0.5)*log(Vm)
        /(2*b_*pow(T,2));
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar soaveRedlichKwong::isobarExpCoef(const scalar rho,const scalar T)  const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar soaveRedlichKwong::isothermalCompressiblity(const scalar rho,const scalar T) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
}


//- Return density [kg/m^3]
inline scalar soaveRedlichKwong::rho(
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
inline scalar soaveRedlichKwong::rho(const scalar p,const scalar T) const
{

    //using perfect gas equation as starting point
    return rho(p,T,p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
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
