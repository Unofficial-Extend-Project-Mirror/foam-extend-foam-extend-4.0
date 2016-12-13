/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    b2_(pr.b2_),
    rhoMin_(pr.rhoMin_),
    rhoMax_(pr.rhoMax_),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0),
    rhostd_(pr.rhostd_)
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar pengRobinson::rhostd() const
{
    return rhostd_;
}


inline scalar pengRobinson::rhoMin() const
{
    return rhoMin_;
}


inline scalar pengRobinson::rhoMax() const
{
    return rhoMax_;
}


inline scalar pengRobinson::Tcrit() const
{
    return Tcrit_;
}


inline void pengRobinson::updateModelCoefficients(const scalar T) const
{
    scalar TTc = sqrt(T/Tcrit_);

    aSave = a0_*(1 + n_*(1 - TTc))*(1 + n_*(1 - TTc));
    daSave = a0_*n_*(n_*TTc - n_ - 1)*TTc/T;
    d2aSave = a0_*n_*(n_ + 1)*TTc/(2*T*T);

    //CL: saving the temperature at which the coefficients are valid
    TSave = T;
}


//CL: Model coefficient a(T)
inline scalar pengRobinson::a(const scalar T) const
{
    //CL: check if a has already been calculated for this temperature
    if(TSave == T)
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
inline scalar pengRobinson::dadT(const scalar T) const
{
    // check if a has already been calculated for this temperature
    if(TSave == T)
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
inline scalar pengRobinson::d2adT2(const scalar T) const
{
    // check if a has already been calculated for this temperature
    if(TSave == T)
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


inline scalar pengRobinson::a0() const
{
    return a0_;
}


inline scalar  pengRobinson::b() const
{
    return b_;
}


inline scalar pengRobinson::n() const
{
    return n_;
}


//returns the pressure for a given density and temperature
inline scalar pengRobinson::p(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    return this->RR()*T/(Vm - b_) - a(T)/(Vm2 + 2*b_*Vm - b2_);
}


//Real deviative dp/dv at constant temperature
//(molar values)
inline scalar pengRobinson::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;

    scalar t1 = bmVm2 - 2*Vm*Vm;
    scalar t12 = t1*t1;

    return (2*a(T)*bmVm2*bpVm - this->RR()*T*t12)/(bmVm2*t12);
}


//Real deviative dp/dT at constant molar volume
//(molar values)
inline scalar pengRobinson::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    return this->RR()/(Vm - b_) - dadT(T)/(Vm2 + 2*b_*Vm - b2_);
}


//Real  deviative dv/dT at constant pressure
//by using implicit differentiation
//(molar values)
inline scalar pengRobinson::dvdT(const scalar rho, const scalar T) const
{
    return -this->dpdT(rho, T)/this->dpdv(rho, T);
}


//Real  deviative dv/dp at constant temperature
//(molar values)
inline scalar pengRobinson::dvdp(const scalar rho, const scalar T)  const
{
    return 1/this->dpdv(rho, T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar pengRobinson::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    scalar root2=pow(2, 0.5);

    return
      - root2*a(T)*log(b_*(1-root2) + Vm)/(4*b_)
      + this->RR()*T*log(Vm - b_)
      + root2*a(T)*log(b_*(root2 + 1) + Vm)/(4*b_);
}


//needed to calculate the entropy
//(molar values)
//needed to calculate the entropy
inline scalar pengRobinson::integral_dpdT_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    scalar root2 = pow(2, 0.5);

    return
      - root2*dadT(T)*log(b_*(1 - root2) + Vm)/(4*b_)
      + this->RR()*log(Vm - b_) + root2*dadT(T)*log(b_*(root2 + 1) + Vm)/(4*b_);
}


//(molar values)
inline scalar pengRobinson::d2pdT2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    return -d2adT2(T)/(Vm2 + 2*b_*Vm - b2_);
}


//(molar values)
inline scalar pengRobinson::d2pdv2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;

    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bmVm3 = bmVm2*bmVm;

    scalar t1 = bmVm2 - 2*Vm2;
    scalar t12 = t1*t1;
    scalar t13 = t12*t1;

    return 2*(a(T)*bmVm3*(5*b_ + 6*b_*Vm + 3*Vm2) - this->RR()*T*t13)/(bmVm*t1);
}


//(molar values)
//using second order implicit differentiation
inline scalar pengRobinson::d2vdT2
(
    const scalar rho,
    const scalar T
) const
{
    scalar dpdT = this->dpdT(rho, T);
    scalar dpdT2 = dpdT*dpdT;
    scalar dpdv = this->dpdv(rho, T);
    scalar dpdv2 = dpdv*dpdv;
    scalar dpdv3 = dpdv2*dpdv;

    return
    -(
          dpdT2*this->d2pdv2(rho, T)
        + dpdv2*this->d2pdT2(rho, T)
        - 2*dpdv*dpdT*this->d2pdvdT(rho, T)
    )
    /dpdv3;
}


//(molar values)
inline scalar pengRobinson::d2pdvdT
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;

    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;

    scalar t1 = bmVm2 - 2*Vm*Vm;
    scalar t12 = t1*t1;
    scalar t14 = t12*t12;

    return (2*dadT(T)*bmVm2*bpVm - this->RR()*t14)/(bmVm2*t12);
}


// the result of this intergal is needed for the nasa based cp polynomial
//(molar values)
inline scalar pengRobinson::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    scalar root2=pow(2, 0.5);

    return root2*d2adT2(T)/(4*b_)
        *(log(b_*(root2 + 1) + Vm) - log(b_*(1 - root2) + Vm));

      //  root2*d2adT2(T)*log(b_*(root2 + 1) + Vm)/(4*b_)
      //- root2*d2adT2(T)*log(b_*(1 - root2)+Vm)/(4*b_);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar pengRobinson::isobarExpCoef
(
    const scalar rho,
    const scalar T
) const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar pengRobinson::isothermalCompressiblity
(
    const scalar rho,
    const scalar T
) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //CL: also possible
    //CL: return -this->dvdp(rho, T)*rho/this->W();
}


//- Return density [kg/m^3]
inline scalar pengRobinson::rho
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
   scalar molarVolumePrevIteration;
   scalar molarVolume;
   label iter = 0;
   label maxIter_ = 400;
   scalar tol_ = 1e-8;
   scalar rho1 = rhoMax_;
   scalar rho2 = rhoMin_;

   molarVolume = this->W()/rho0;

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
                   /(this->dpdv((this->W()/molarVolumePrevIteration), T))
               )/pow(2,i);

           i++;
           if (i>8)
           {
               //CL: using bisection methode as backup,
               //CL: solution must be between rhoMin_ to rhoMax
               for(i=0; i<200; i++)
               {
                   scalar f1 = this->p(rho1, T) - p;
                   scalar f2 = this->p(rho2, T) - p;
                   scalar rho3 = (rho1 + rho2)/2;
                   scalar f3 = this->p(rho3, T) - p;

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
           mag(this->p((this->W()/molarVolume), T) - p)
         > mag(this->p((this->W()/molarVolumePrevIteration), T) - p)
       );

       if (iter++ > maxIter_)
       {
           FatalErrorIn
           (
               "inline scalar pengRobinson::rho"
               "(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

   return this->W()/molarVolume;
}


//- Return density [kg/m^3]on
inline scalar pengRobinson::rho(const scalar p, const scalar T) const
{
    //CL: using perfect gas equation as starting point
    return rho(p, T, p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar pengRobinson::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho, T)*rho*rho/this->W();
}


//- Return compression factor []
inline scalar pengRobinson::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    return p/(this->R()*T*this->rho(p, T, rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void pengRobinson::operator+=(const pengRobinson& pr)
{
    specie::operator+=(pr);
}


inline void pengRobinson::operator*=(const scalar s)
{
    specie::operator*=(s);
}


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
