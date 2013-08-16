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

#include "aungierRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline aungierRedlichKwong::aungierRedlichKwong
(
    const specie& sp
)
:
    specie(sp)
{}


// Construct as named copy
inline aungierRedlichKwong::aungierRedlichKwong(const word& name, const aungierRedlichKwong& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
    azentricFactor_(pg.azentricFactor_),
    rhocrit_(pg.rhocrit_),    
    a0_(pg.a0_),
    b_(pg.b_),
    c_(pg.c_),
    n_(pg.n_),
    rhostd_(pg.rhostd_),
    b2_(pg.b2_),
    b3_(pg.b3_),
    b4_(pg.b4_),
    b5_(pg.b5_),
    c2_(pg.c2_)
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

inline void aungierRedlichKwong::updateModelCoefficients(const scalar T)const
{
    aSave=a0_*pow(T/Tcrit_,-n());
    daSave=-a0_*n()*pow(T/Tcrit_,-n())/T;
    d2aSave=a0_*(n() *n()+n())/(T*T)*pow(T/Tcrit_,-n());

    //CL: saving the temperature at which the coefficients are valid
    TSave=T;
}


//CL: Model coefficient a(T)
inline scalar aungierRedlichKwong::a(const scalar T)const
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
inline scalar aungierRedlichKwong::dadT(const scalar T)const
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
inline scalar aungierRedlichKwong::d2adT2(const scalar T)const
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

//Aungier Redlich Kwong factors
inline scalar aungierRedlichKwong::a0()const
{
    return a0_;
}


inline scalar aungierRedlichKwong::b()const
{
    return b_;
}

inline scalar aungierRedlichKwong::c()const
{
    return c_;
}

inline scalar aungierRedlichKwong::n()const
{
    return n_;
}

//CL: pow of constants (b(), c()) used in the code e.g. b2_=b*b;
inline scalar aungierRedlichKwong::b2()const
{
    return b2_;
}

inline scalar aungierRedlichKwong::b3()const
{
    return b3_;
}

inline scalar aungierRedlichKwong::b4()const
{
    return b4_;
}

inline scalar aungierRedlichKwong::b5()const
{
    return b5_;
}

inline scalar aungierRedlichKwong::c2()const
{
    return c2_;
}


//returns the pressure for a given density and temperature    
inline scalar aungierRedlichKwong::p(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return this->RR*T/(Vm-b()+c())
               -a(T)/(Vm*(Vm+b()));
}


//Real deviative dp/dv at constant temperature 
//(molar values)
inline scalar aungierRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
   scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    return  (a(T)*(b3()-2*b2()*c()+b()*(c()+Vm)*(c()-3*Vm)+2*Vm*pow(c()+Vm,2))
                -this->RR*T*Vm2*(b2()+2*b()*Vm+Vm2))
                    /(Vm2*pow(b()-c()-Vm,2)*pow(b()+Vm,2));
}


//Real deviative dp/dT at constant molar volume 
//(molar values)
inline scalar aungierRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return this->RR/(Vm-b()+c())-dadT(T)/(Vm*(Vm+b()));
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
    return this->RR*T*log(-b()+c()+Vm)+a(T)*log(b()+Vm)/b()-a(T)*log(Vm)/b();
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
    return this->RR*log(-b()+c()+Vm)+dadT(T)*log(b()+Vm)/b()-dadT(T)*log(Vm)/b();
}


//(molar values)
inline scalar aungierRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    return -d2adT2(T)/(Vm*(Vm+b()));
}


//(molar values)
inline scalar aungierRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm*Vm*Vm;

    return -2*(a(T)*(b5()-3*b4()*c()+3*b3()*(c2()-c()*Vm-Vm2)
               -b2()*(c()+Vm)*(c2()-7*c()*Vm+Vm2)+3*b()*Vm*pow(c()+Vm,2)*(2*Vm-c())
                   -3*Vm2*pow(c()+Vm,3))+this->RR*T*Vm3*(b3()+3*b2()*Vm
                        +3*b()*Vm2+Vm3))
                            /(Vm3*pow(b()-c()-Vm,3)*pow(b()+Vm,3));
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
    scalar Vm2 = Vm*Vm;

    return (dadT(T)*(b3()-2*b2()*c()+b()*(c()+Vm)*(c()-3*Vm)+2*Vm*pow(c()+Vm,2))
               -this->RR*Vm2*(b2()+2*b()*Vm+Vm2))
                   /(Vm2*pow(b()-c()-Vm,2)*pow(b()+Vm,2));
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
    return d2adT2(T)*log(b()+Vm)/b()-d2adT2(T)*log(Vm)/b();
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
               "inline scalar aungierRedlichKwong::rho(const scalar p, const scalar T, const scalar rho0) const "
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
    return p/(this->R()*T*this->rho(p,T,rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void aungierRedlichKwong::operator+=(const aungierRedlichKwong& pr)
{
    specie::operator+=(pr);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline aungierRedlichKwong operator+
(
    const aungierRedlichKwong& pr1,
    const aungierRedlichKwong& pr2
)
{
    return aungierRedlichKwong
    (
        static_cast<const specie&>(pr1)
      + static_cast<const specie&>(pr2)
    );
}

inline aungierRedlichKwong operator*
(
    const scalar s,
    const aungierRedlichKwong& pr
)
{
    return aungierRedlichKwong(s*static_cast<const specie&>(pr));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
