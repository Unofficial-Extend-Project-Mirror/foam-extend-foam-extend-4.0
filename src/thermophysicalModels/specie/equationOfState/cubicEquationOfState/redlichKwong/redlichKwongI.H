/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "redlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline redlichKwong::redlichKwong
(
    const specie& sp
)
:
    specie(sp)
{}


// Construct as named copy
inline redlichKwong::redlichKwong(const word& name, const redlichKwong& rk)
:
    specie(name, rk),
    pcrit_(rk.pcrit_),
    Tcrit_(rk.Tcrit_),
    a_(rk.a_),
    b_(rk.b_),
    rhostd_(rk.rhostd_),
    b2_(rk.b2_),
    b3_(rk.b3_),
    b5_(rk.b5_)
{}


// Construct and return a clone
inline autoPtr<redlichKwong> redlichKwong::clone() const
{
    return autoPtr<redlichKwong>(new redlichKwong(*this));
}


// Selector from Istream
inline autoPtr<redlichKwong> redlichKwong::New(Istream& is)
{
    return autoPtr<redlichKwong>(new redlichKwong(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar redlichKwong::rhostd()const
{
    return rhostd_;
}


inline scalar redlichKwong::a()const
{
    return a_;
}


inline scalar redlichKwong::b()const
{
    return b_;
}

inline scalar redlichKwong::b2()const
{
    return b2_;
}

inline scalar redlichKwong::b3()const
{
    return b3_;
}

inline scalar redlichKwong::b5()const
{
    return b5_;
}



//returns the pressure for a given density and temperature    
inline scalar redlichKwong::p(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;    
    return this->RR*T/(Vm - b_) - a_/(sqrt(T)*Vm*(Vm + b_));
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar redlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    return (a_*(b3() - 3*b_*Vm2 + 2*pow(Vm,3))
    - this->RR*pow(T,1.5)*Vm2*(b2() + 2*b_*Vm + Vm2))
        /(sqrt(T)*Vm2*pow((b_ + Vm),2)*pow( (b_ - Vm),2));
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar redlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return 0.5*a_/(pow(T,1.5)*Vm*(b_ + Vm))-this->RR/(b_ - Vm);
}


//Real  deviative dv/dT at constant pressure 
//using implicit differentiation
//(molar values)
inline scalar redlichKwong::dvdT(const scalar rho, const scalar T) const
{
    return -this->dpdT(rho,T)/this->dpdv(rho,T);
}


//Real  deviative dv/dp at constant temperature 
//(molar values)
inline scalar redlichKwong::dvdp(const scalar rho, const scalar T)  const 
{
    return 1/this->dpdv(rho,T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar redlichKwong::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    return this->RR*T*log(Vm - b_)
        + (a_*log(b_ + Vm))/(b_*sqrt(T))
        - (a_*log(Vm))/(b_*sqrt(T));
}


//needed to calculate the entropy
//(molar values)
inline scalar redlichKwong::integral_dpdT_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    return this->RR*log(Vm - b_)
        -(a_*log(b_ + Vm))/(2*b_*pow(T,1.5))
        +(a_*log(Vm))/(2*b_*pow(T,1.5));
}


//(molar values)
inline scalar redlichKwong::d2pdT2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -0.75*a_/(pow(T,2.5)*Vm*(b_ + Vm));
}


//(molar values)
inline scalar redlichKwong::d2pdv2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;

    return
    (
        2*(
            a_*(
                b5()-3*b3()*Vm2
                - b2()*Vm3
                + 6*b_*pow(Vm,4)-3*pow(Vm,5)
            )
            + this->RR*pow(T,1.5)*Vm3*(b3()
            + 3*b2()*Vm
            + 3*b_*Vm2+Vm3)
        )
        /(sqrt(T)*Vm3*pow((b_ + Vm),3)*pow(Vm-b_,3))
    );
}


//(molar values)
//using second Order implicit differentiation 
inline scalar redlichKwong::d2vdT2(const scalar rho, const scalar T) const
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
inline scalar redlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    
    return
    -(0.5*(
        a_*(b3() - 3*b_*Vm2 + 2*pow(Vm,3))
      + 2*this->RR*pow(T,1.5)*Vm2*(b2() + 2*b_*Vm + Vm2)
    ))
    /(pow(T,1.5)*Vm2*pow(b_ + Vm,2)*pow(b_ - Vm,2));
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar redlichKwong::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
) const 
{
    scalar Vm = this->W()/rho;
    return 0.75*a_*log(b_  +  Vm)/(pow(T,2.5)*b_)
        - 0.75*a_*log(Vm)/(pow(T,2.5)*b_);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar redlichKwong::isobarExpCoef
(
    const scalar rho,
    const scalar T
)  const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar redlichKwong::isothermalCompressiblity
(
    const scalar rho,
    const scalar T
) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //also possible : return -this->dvdp(rho,T)*rho/this->W();
}


//- Return density [kg/m^3]
inline scalar redlichKwong::rho
(
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
           //CL: modified Newton solver
           molarVolume=molarVolumePrevIteration
               -(
                   (this->p((this->W()/molarVolumePrevIteration),T) - p)
                   /(this->dpdv((this->W()/molarVolumePrevIteration),T))
               )/pow(2,i);

           i++;
           if (i>8) 
           {
               //CL: using bisection methode as backup,
               //CL: solution must be between rhoMin_ to rhoMax
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
inline scalar redlichKwong::rho(const scalar p, const scalar T) const
{
    // using perfect gas equation as starting point
    return rho(p,T,p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar redlichKwong::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho,T)*pow(rho,2)/this->W();
}


//- Return compression factor []
inline scalar redlichKwong::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    return p/(this->R()*T*this->rho(p,T,rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void redlichKwong::operator+=(const redlichKwong& rk)
{
    specie::operator+=(rk);
}

/*
inline void redlichKwong::operator-=(const redlichKwong& rk)
{
    specie::operator-=(rk);
}


inline void redlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}

*/
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


inline redlichKwong operator+
(
    const redlichKwong& rk1,
    const redlichKwong& rk2
)
{
    return redlichKwong
    (
        static_cast<const specie&>(rk1)
      + static_cast<const specie&>(rk2)
    );
}

inline redlichKwong operator*
(
    const scalar s,
    const redlichKwong& rk
)
{
    return redlichKwong(s*static_cast<const specie&>(rk));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
