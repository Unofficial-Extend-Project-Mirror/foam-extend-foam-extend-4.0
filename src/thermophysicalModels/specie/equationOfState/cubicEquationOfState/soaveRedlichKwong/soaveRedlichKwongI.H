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
inline soaveRedlichKwong::soaveRedlichKwong(const word& name, const soaveRedlichKwong& srk)
:
    specie(name, srk),
    pcrit_(srk.pcrit_),
    Tcrit_(srk.Tcrit_),
    azentricFactor_(srk.azentricFactor_),
    a0_(srk.a0_),
    b_(srk.b_),
    n_(srk.n_),
    b2_(srk.b2_),
    rhoMin_(srk.rhoMin_),
    rhoMax_(srk.rhoMax_),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0),
    rhostd_(srk.rhostd_)
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

inline scalar soaveRedlichKwong::rhostd() const
{
    return rhostd_;
}


inline scalar soaveRedlichKwong::rhoMin() const
{
    return rhoMin_;
}


inline scalar soaveRedlichKwong::rhoMax() const
{
    return rhoMax_;
}


inline scalar soaveRedlichKwong::Tcrit() const
{
    return Tcrit_;
}


inline void soaveRedlichKwong::updateModelCoefficients(const scalar T) const
{
    scalar TTc = pow(T/Tcrit_, 0.5);

    aSave = a0_*(1 + n_*(1 - TTc))*(1 + n_*(1 - TTc));
    daSave = a0_*n_*(n_*TTc - n_ - 1)*TTc/T;
    d2aSave = a0_*n_*(n_ + 1)*TTc/(2*T*T);

    //CL: saving the temperature at which the coefficients are valid
    TSave = T;
}


//CL: Model coefficient a(T)
inline scalar soaveRedlichKwong::a(const scalar T) const
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
inline scalar soaveRedlichKwong::dadT(const scalar T) const
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
inline scalar soaveRedlichKwong::d2adT2(const scalar T) const
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


inline scalar soaveRedlichKwong::a0() const
{
    return a0_;
}


inline scalar soaveRedlichKwong::b() const
{
    return b_;
}


inline scalar soaveRedlichKwong::n() const
{
    return n_;
}


//returns the pressure for a given density and temperature
inline scalar soaveRedlichKwong::p(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    return this->RR()*T/(Vm - b_) - a(T)/(Vm*(Vm + b_));
}


//Real deviative dp/dv at constant temperature
//(molar values)
inline scalar soaveRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;

    return(a(T)*bmVm2*(b_ + 2*Vm) - this->RR()*T*Vm2*bpVm2)/(Vm2*bpVm2*bmVm2);
}


//Real deviative dp/dT at constant molar volume
//(molar values)
inline scalar soaveRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR()/(Vm - b_) - dadT(T)/(Vm*(Vm + b_));
}


//Real  deviative dv/dT at constant pressure
//using implicit differentiation
// (molar values)
inline scalar soaveRedlichKwong::dvdT(const scalar rho, const scalar T) const
{
    return -this->dpdT(rho, T)/this->dpdv(rho, T);
}


//Real  deviative dv/dp at constant temperature
//(molar values)
inline scalar soaveRedlichKwong::dvdp(const scalar rho, const scalar T) const
{
    return 1/this->dpdv(rho, T);
}


//needed to calculate the internal energy
//(molar values)
inline scalar soaveRedlichKwong::integral_p_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;

    return this->RR()*T*log(Vm - b_) + a(T)/b_*log((b_ + Vm)/Vm);
}


//needed to calculate the entropy
//(molar values)
//needed to calculate the entropy
inline scalar soaveRedlichKwong::integral_dpdT_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;

    return this->RR()*log(Vm - b_) + dadT(T)/b_*log((b_ + Vm)/Vm);
}


//(molar values)
inline scalar soaveRedlichKwong::d2pdT2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -d2adT2(T)/(Vm*(Vm + b_));
}


//(molar values)
inline scalar soaveRedlichKwong::d2pdv2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm2*Vm;
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bmVm3 = bmVm2*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;
    scalar bpVm3 = bpVm2*bpVm;

    return 2*(a(T)*bmVm3*(b2_ + 3*b_*Vm + 3*Vm2) + this->RR()*T*Vm3*bpVm3)
        /(Vm3*bpVm3*bmVm3);
}


//(molar values)
//using second order implicit differentiation
inline scalar soaveRedlichKwong::d2vdT2
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
inline scalar soaveRedlichKwong::d2pdvdT
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;

    return (dadT(T)*bmVm2*(b_ + 2*Vm) - this->RR()*Vm2*bmVm2)
        /(Vm2*bpVm2*bmVm2);
}


// the result of this intergal is needed for the nasa based cp polynomial
//(molar values)
inline scalar soaveRedlichKwong::integral_d2pdT2_dv
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;

    return d2adT2(T)/b_*log((b_ +  Vm)/Vm);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar soaveRedlichKwong::isobarExpCoef
(
    const scalar rho,
    const scalar T
) const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar soaveRedlichKwong::isothermalCompressiblity
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
inline scalar soaveRedlichKwong::rho
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
               "inline scalar soaveRedlichKwong::rho"
               "(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

   return this->W()/molarVolume;
}


//- Return density [kg/m^3]
inline scalar soaveRedlichKwong::rho(const scalar p, const scalar T) const
{
    //CL: using perfect gas equation as starting point
    return rho(p, T, p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar soaveRedlichKwong::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho, T)*rho*rho/this->W();
}


//- Return compression factor []
inline scalar soaveRedlichKwong::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    return p/(this->R()*T*this->rho(p, T, rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void soaveRedlichKwong::operator+=(const soaveRedlichKwong& srk)
{
    specie::operator+=(srk);
}


inline void soaveRedlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}


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
