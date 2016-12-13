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
    specie(sp),
    TSave(0)
{}


// Construct as named copy
inline aungierRedlichKwong::aungierRedlichKwong(const word& name, const aungierRedlichKwong& pg)
:
    specie(name, pg),
    pcrit_(pg.pcrit_),
    Tcrit_(pg.Tcrit_),
    rhocrit_(pg.rhocrit_),
    azentricFactor_(pg.azentricFactor_),
    a0_(pg.a0_),
    b_(pg.b_),
    c_(pg.c_),
    n_(pg.n_),
    b2_(pg.b2_),
    b3_(pg.b3_),
    b4_(pg.b4_),
    b5_(pg.b5_),
    c2_(pg.c2_),
    rhoMin_(pg.rhoMin_),
    rhoMax_(pg.rhoMax_),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0),
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar aungierRedlichKwong::rhostd() const
{
    return rhostd_;
}


inline scalar aungierRedlichKwong::rhoMin() const
{
    return rhoMin_;
}


inline scalar aungierRedlichKwong::rhoMax() const
{
    return rhoMax_;
}


inline scalar aungierRedlichKwong::Tcrit() const
{
    return Tcrit_;
}


inline scalar aungierRedlichKwong::pcrit() const
{
    return pcrit_;
}


inline scalar aungierRedlichKwong::rhocrit() const
{
    return rhocrit_;
}


inline scalar aungierRedlichKwong::azentricFactor() const
{
    return azentricFactor_;
}


inline void aungierRedlichKwong::updateModelCoefficients(const scalar T)const
{
    scalar TTc = pow(T/Tcrit_, -n_);

    aSave = a0_*TTc;
    daSave = -a0_*n_*TTc/T;
    d2aSave = a0_*(n_*n_ + n_)/(T*T)*TTc;

    //CL: saving the temperature at which the coefficients are valid
    TSave = T;
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
inline scalar aungierRedlichKwong::a0() const
{
    return a0_;
}


inline scalar aungierRedlichKwong::b() const
{
    return b_;
}


inline scalar aungierRedlichKwong::c() const
{
    return c_;
}

inline scalar aungierRedlichKwong::n() const
{
    return n_;
}


//returns the pressure for a given density and temperature
inline scalar aungierRedlichKwong::p(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR()*T/(Vm - b_ + c_) - a(T)/(Vm*(Vm + b_));
}


//Real deviative dp/dv at constant temperature
//(molar values)
inline scalar aungierRedlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;
    scalar cpVm = c_ + Vm;

    return
    (
        a(T)*(b3_ - 2*b2_*c_ + b_*cpVm*(c_ - 3*Vm) + 2*Vm*cpVm*cpVm)
      - this->RR()*T*Vm2*bpVm2
    )
   /(Vm2*pow(b_ - c_ - Vm, 2)*bpVm2);
}


//Real deviative dp/dT at constant molar volume
//(molar values)
inline scalar aungierRedlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR()/(Vm - b_ + c_)-dadT(T)/(Vm*(Vm + b_));
}


//Real  deviative dv/dT at constant pressure
//using implicit differentiation
//(molar values)
inline scalar aungierRedlichKwong::dvdT(const scalar rho, const scalar T) const
{
    return -this->dpdT(rho, T)/this->dpdv(rho, T);
}


//Real  deviative dv/dp at constant temperature
//(molar values)
inline scalar aungierRedlichKwong::dvdp(const scalar rho, const scalar T) const
{
    return 1/this->dpdv(rho, T);
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

    return this->RR()*T*log(Vm - b_ + c_) + a(T)/b_*log((b_ + Vm)/Vm);
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

    return this->RR()*log(Vm - b_ + c_) + dadT(T)/b_*log((b_ + Vm)/Vm);
}


//(molar values)
inline scalar aungierRedlichKwong::d2pdT2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -d2adT2(T)/(Vm*(b_ + Vm));
}


//(molar values)
inline scalar aungierRedlichKwong::d2pdv2(const scalar rho,const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm2*Vm;
    scalar cpVm = c_ + Vm;
    scalar bpVm3 = pow(b_ + Vm, 3);

    return -2*
    (
        a(T)*
        (
            b5_ - 3*b4_*c_ + 3*b3_*pow(c_ - Vm, 2)
          + cpVm*
            (
              - b2_*(c2_ - 7*c_*Vm + Vm2)
              + cpVm*
                (
                    3*b_*Vm*(2*Vm - c_)
                  - 3*Vm2*cpVm
                )
            )
        )
      + this->RR()*T*Vm3*bpVm3
    )
   /(Vm3*pow(b_ - c_ - Vm, 3)*bpVm3);
}


//(molar values)
//using second order implicit differentiation
inline scalar aungierRedlichKwong::d2vdT2
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
inline scalar aungierRedlichKwong::d2pdvdT
(
    const scalar rho,
    const scalar T
) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar cpVm = c_ + Vm;
    scalar cpVm2 = cpVm*cpVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;

    return
    (
        dadT(T)*(b3_ - 2*b2_*c_ + b_*cpVm*(c_ - 3*Vm) + 2*Vm*cpVm2)
      - this->RR()*Vm2*bpVm2
    )
    /(Vm2*pow(b_ - c_ - Vm, 2)*bpVm2);
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

    return d2adT2(T)/b_*log((b_ + Vm)/Vm);
}


//Isobar expansion Coefficent beta = 1/v (dv/dt) at constant p
//(molar values)
inline scalar aungierRedlichKwong::isobarExpCoef
(
    const scalar rho,
    const scalar T
) const
{
    return this->dvdT(rho, T)*rho/this->W();
}


//isothemal compressiblity kappa (not Thermal conductivity)
//(molar values)
inline scalar aungierRedlichKwong::isothermalCompressiblity
(
    const scalar rho,
    const scalar T
) const
{
    return this->isobarExpCoef(rho, T)/this->dpdT(rho, T);
    //also possible : return -this->dvdp(rho,T)*rho/this->W();
}


//- Return density [kg/m^3]
inline scalar aungierRedlichKwong::rho
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
               "inline scalar aungierRedlichKwong::rho"
               "(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

   return this->W()/molarVolume;
}


//- Return density [kg/m^3]
inline scalar aungierRedlichKwong::rho(const scalar p, const scalar T) const
{
    //CL: using perfect gas equation as starting point
    return rho(p, T, p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar aungierRedlichKwong::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho, T)*rho*rho/this->W();
}


//- Return compression factor []
inline scalar aungierRedlichKwong::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    return p/(this->R()*T*this->rho(p, T, rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void aungierRedlichKwong::operator+=(const aungierRedlichKwong& ark)
{
    specie::operator+=(ark);
}


inline void aungierRedlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline aungierRedlichKwong operator+
(
    const aungierRedlichKwong& ark1,
    const aungierRedlichKwong& ark2
)
{
    return aungierRedlichKwong
    (
        static_cast<const specie&>(ark1)
      + static_cast<const specie&>(ark2)
    );
}


inline aungierRedlichKwong operator*
(
    const scalar s,
    const aungierRedlichKwong& ark
)
{
    return aungierRedlichKwong(s*static_cast<const specie&>(ark));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
