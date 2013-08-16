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

#include "pengRobinson.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline pengRobinson::pengRobinson
(
    const specie& sp    
)
:
    specie(sp),
    TSave(0) 
{}


// Construct as named copy
inline pengRobinson::pengRobinson(const word& name, const pengRobinson& pr)
:
    specie(name, pr),
    pcrit_(pr.pcrit_),
    Tcrit_(pr.Tcrit_),
    azentricFactor_(pr.azentricFactor_),
    a0_(pr.a0_),
    b_(pr.b_),
    n_(pr.n_),
    rhostd_(pr.rhostd_),
    TSave(0),
    b2_(pr.b2_),
    b3_(pr.b3_),
    b4_(pr.b4_),
    b5_(pr.b5_)
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

inline void pengRobinson::updateModelCoefficients(const scalar T)const
{
    aSave=a0_*pow(1+n_*(1-pow(T/Tcrit_,0.5)),2);
    daSave=a0_*n_*(n_*sqrt(T/Tcrit_)-n_-1)*sqrt(T/Tcrit_)/T;
    d2aSave=a0_*n_*(n_+1)*sqrt(T/Tcrit_)/(2*pow(T,2));

    //CL: saving the temperature at which the coefficients are valid
    TSave=T;
}


inline scalar pengRobinson::rhostd()const
{
    return rhostd_;
}


//CL: Model coefficient a(T)
inline scalar pengRobinson::a(const scalar T)const
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
inline scalar pengRobinson::dadT(const scalar T)const
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
inline scalar pengRobinson::d2adT2(const scalar T)const
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


inline scalar pengRobinson::a0()const
{
    return a0_;
}


inline scalar pengRobinson::b()const
{
    return b_;
}


inline scalar pengRobinson::n()const
{
    return n_;
}

//CL: pow of constant b() used in the code e.g. b2_=b*b;
inline scalar pengRobinson::b2()const
{
    return b2_;
}

inline scalar pengRobinson::b3()const
{
    return b3_;
}

inline scalar pengRobinson::b4()const
{
    return b4_;
}

inline scalar pengRobinson::b5()const
{
    return b5_;
}

//returns the pressure for a given density and temperature    
inline scalar pengRobinson::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
     
    return this->RR*T/(Vm-b())-a(T)/(pow(Vm,2)+2*b()*Vm-b2());
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar pengRobinson::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;

    return(
              2*a(T)*
              (
                  b3()-b2()*Vm-b()*Vm2+Vm3
              )
              -this->RR*T*
              (
                  b4()-4*b3()*Vm+2*b2()*Vm2
                      +4*b()*Vm3+pow(Vm,4)
              )
          )
          /(pow(b()-Vm,2)*pow(b2()-2*b()*Vm-Vm2,2));
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar pengRobinson::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR/(Vm-b())-dadT(T)/(pow(Vm,2)+2*b()*Vm-b2());
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

    return -pow(2,0.5)*a(T)*log(b()*(1-pow(2,0.5))+Vm)/(4*b())+this->RR*T*log(Vm-b())
               +pow(2,0.5)*a(T)*log(b()*(pow(2,0.5)+1)+Vm)/(4*b());
}


//(molar values)
//needed to calculate the entropy
inline scalar pengRobinson::integral_dpdT_dv(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -pow(2,0.5)*dadT(T)*log(b()*(1-pow(2,0.5))+Vm)/(4*b())
               +this->RR*log(Vm-b())+pow(2,0.5)*dadT(T)*log(b()*(pow(2,0.5)+1)+Vm)/(4*b());
}


//(molar values)
inline scalar pengRobinson::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -d2adT2(T)/(pow(Vm,2)+2*b()*Vm-b2());
}


//(molar values)
inline scalar pengRobinson::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;
    scalar Vm4 = Vm*Vm*Vm*Vm;
    scalar Vm5 = Vm*Vm*Vm*Vm*Vm;

    return 2*
           (
               a(T)*
               (
                   5*b5()-9*b4()*Vm+4*b2()*Vm3+3*b()*Vm4-3*Vm5
               )
               -this->RR*T*
               (
                   pow(b(),6)-6*b5()*Vm+9*b4()*Vm2+4*b3()*Vm3
                       -9*b2()*Vm4-6*b()*Vm5-pow(Vm,6)
               )
           )
           /(pow(b()-Vm,3)*pow(b2()-2*b()*Vm-Vm2,3));
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
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;
  
    return(
              2*dadT(T)*
              (
                  b3()-b2()*Vm-b()*Vm2+Vm3
              )
              -this->RR*
              (
                  b4()-4*b3()*Vm+2*b2()*Vm2
                      +4*b()*Vm3+pow(Vm,4)
              )
          )
          /(pow(b()-Vm,2)*pow(b2()-2*b()*Vm-Vm2,2));
}


// the result of this intergal is needed for the nasa based cp polynomial 
//(molar values)
inline scalar pengRobinson::integral_d2pdT2_dv(const scalar rho,const scalar T) const 
{
    scalar Vm = this->W()/rho;

    return pow(2,0.5)*d2adT2(T)*log(b()*(pow(2,0.5)+1)+Vm)/(4*b())-pow(2,0.5)*d2adT2(T)*log(b()*(1-pow(2,0.5))+Vm)/(4*b());
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
    //CL: also possible 
    //CL: return -this->dvdp(rho,T)*rho/this->W();
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
               "inline scalar pengRobinson::rho(const scalar p, const scalar T, const scalar rho0) const "
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
    return p/(this->R()*T*this->rho(p,T,rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void pengRobinson::operator+=(const pengRobinson& pr)
{
    specie::operator+=(pr);
}

/*
inline void pengRobinson::operator-=(const pengRobinson& pr)
{
    specie::operator-=(pr);
}

inline void pengRobinson::operator*=(const scalar s)
{
    specie::operator*=(s);
}

*/
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline pengRobinson operator+
(
    const pengRobinson& pr1,
    const pengRobinson& pr2
)
{
    return pengRobinson
    (
        static_cast<const specie&>(pr1)
      + static_cast<const specie&>(pr2)
    );
}

inline pengRobinson operator*
(
    const scalar s,
    const pengRobinson& pr
)
{
    return pengRobinson(s*static_cast<const specie&>(pr));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
