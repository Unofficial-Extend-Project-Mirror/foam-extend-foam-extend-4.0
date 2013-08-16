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

#include "soaveRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline soaveRedlichKwong::soaveRedlichKwong
(
    const specie& sp   
)
:
    specie(sp),
    TSave(0)     
{}


// Construct as named copy
inline soaveRedlichKwong::soaveRedlichKwong(const word& name,const  soaveRedlichKwong& srk)
:
    specie(name, srk),
    pcrit_(srk.pcrit_),
    Tcrit_(srk.Tcrit_),
    azentricFactor_(srk.azentricFactor_),
    a0_(srk.a0_),
    b_(srk.b_),
    n_(srk.n_),
    rhostd_(srk.rhostd_),
    TSave(0),
    b2_(srk.b2_),
    b3_(srk.b3_),
    b5_(srk.b5_) 
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

inline void soaveRedlichKwong::updateModelCoefficients(const scalar T)const
{
    aSave=a0_*pow(1+n_*(1-pow(T/Tcrit_,0.5)),2);
    daSave=a0_*n_*(n_*sqrt(T/Tcrit_)-n_-1)*sqrt(T/Tcrit_)/T;
    d2aSave=a0_*n_*(n_+1)*sqrt(T/Tcrit_)/(2*pow(T,2));

    //CL: saving the temperature at which the coefficients are valid
    TSave=T;
}


inline scalar soaveRedlichKwong::rhostd()const
{
    return rhostd_;
}


//CL: Model coefficient a(T)
inline scalar soaveRedlichKwong::a(const scalar T)const
{
    //CL: check if a has already been calculated for this temperature
    if(TSave==T)
    { 
        return aSave;  
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return aSave;
    }     
}


//CL: temperature deriviative of model coefficient a(T)
inline scalar soaveRedlichKwong::dadT(const scalar T)const
{
    // check if a has already been calculated for this temperature
    if(TSave==T)
    {
        return daSave;  
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return daSave;
    }    
}


//CL: second order temperature deriviative of model coefficient a(T)
inline scalar soaveRedlichKwong::d2adT2(const scalar T)const
{
    // check if a has already been calculated for this temperature
    if(TSave==T)
    {
        return d2aSave;  
    }
    //CL: If not, recalculate a(T), dadT(T) and d2adT2(T)
    else
    {
        updateModelCoefficients(T);
        return d2aSave;
    }    
}


inline scalar soaveRedlichKwong::a0()const
{
    return a0_;
}


inline scalar soaveRedlichKwong::b()const
{
    return b_;
}


inline scalar soaveRedlichKwong::n()const
{
    return n_;
}

//CL: pow of constant b() used in the code e.g. b2_=b*b;
inline scalar soaveRedlichKwong::b2()const
{
    return b2_;
}

inline scalar soaveRedlichKwong::b3()const
{
    return b3_;
}

inline scalar soaveRedlichKwong::b5()const
{
    return b5_;
}


//returns the pressure for a given density and temperature    
inline scalar soaveRedlichKwong::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return 
    (
        this->RR*T/(Vm-b_)
        -a(T)/(Vm*(Vm+b_))
    );
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar soaveRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    return 
    (
        a(T)*(b3()-3*b_*Vm2+2*pow(Vm,3))
        -this->RR*T*Vm2*(b2()+2*b_*Vm+Vm2)
    )
    /(Vm2*pow(b_+Vm,2)*pow(b_-Vm,2));   
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar soaveRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
 
    return 
    (
        this->RR/(Vm-b_)
        -dadT(T)/(Vm*(Vm+b_))
    );
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
 
    return this->RR*T*log(Vm-b_)+a(T)*log(b_+Vm)/b_-a(T)*log(Vm)/b_;
}


//needed to calculate the entropy
//(molar values)
inline scalar soaveRedlichKwong::integral_dpdT_dv(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR*log(Vm-b_)+dadT(T)*log(b_+Vm)/b_-dadT(T)*log(Vm)/b_;
}


//(molar values)
inline scalar soaveRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -d2adT2(T)/(Vm*(Vm+b_));
}


//(molar values) 
inline scalar soaveRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;

    return 
    2*
    (
        a(T)*
        (
            b5()-3*b3()*Vm2-b2()*Vm3+6*b_*pow(Vm,4)-3*pow(Vm,5)
        )
        +this->RR*T*Vm3*
        (
            b3()+3*b2()*Vm+3*b_*Vm2+Vm3
        )
    )
    /(Vm3*pow(b_+Vm,3)*pow(Vm-b_,3));   
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
    scalar Vm2 = Vm*Vm;

    return 
    (
        dadT(T)*(b3()-3*b_*Vm2+2*pow(Vm,3))
        -this->RR*Vm2*(b2()+2*b_*Vm+Vm2)
    )
    /(Vm2*pow(b_+Vm,2)*pow(b_-Vm,2)); 
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar soaveRedlichKwong::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{    
    scalar Vm = this->W()/rho;

    return d2adT2(T)*log(b_+Vm)/b_-d2adT2(T)*log(Vm)/b_;
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
               //CL: using bisection methode as backup,
               //CL: solution must be between rho=0.001 to rho=1500;
               //CL: if not, change rhoMax_ and rhoMin_
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
               "inline scalar soaveRedlichKwong::rho(const scalar p, const scalar T, const scalar rho0) const "
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
    //CL: using perfect gas equation as starting point
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
    return p/(this->R()*T*this->rho(p,T,rho0));
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void soaveRedlichKwong::operator+=(const soaveRedlichKwong& srk)
{
    specie::operator+=(srk);
}

/*
inline void soaveRedlichKwong::operator-=(const soaveRedlichKwong& srk)
{
    specie::operator-=(srk);
}

inline void soaveRedlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}

*/
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline soaveRedlichKwong operator+
(
    const soaveRedlichKwong& srk1,
    const soaveRedlichKwong& srk2
)
{
    return soaveRedlichKwong
    (
        static_cast<const specie&>(srk1)
      + static_cast<const specie&>(srk2)
    );
}

inline soaveRedlichKwong operator*
(
    const scalar s,
    const soaveRedlichKwong& srk
)
{
    return soaveRedlichKwong(s*static_cast<const specie&>(srk));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
