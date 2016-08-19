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
    b2_(rk.b2_),
    b3_(rk.b3_),
    rhoMin_(rk.rhoMin_),
    rhoMax_(rk.rhoMax_),
    rhostd_(rk.rhostd_)
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

inline scalar redlichKwong::rhostd() const
{
    return rhostd_;
}


inline scalar redlichKwong::rhoMin() const
{
    return rhoMin_;
}


inline scalar redlichKwong::rhoMax() const
{
    return rhoMax_;
}


inline scalar redlichKwong::Tcrit() const
{
    return Tcrit_;
}


inline scalar redlichKwong::a() const
{
    return a_;
}


inline scalar redlichKwong::b() const
{
    return b_;
}


//returns the pressure for a given density and temperature
inline scalar redlichKwong::p(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return this->RR()*T/(Vm - b_) - a_/(sqrt(T)*Vm*(Vm + b_));
}


//Real deviative dp/dv at constant temperature
//(molar values)
inline scalar redlichKwong::dpdv(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm2*Vm;
    scalar T05 = pow(T, 0.5);
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;

    return (a_*(b3_ - 3*b_*Vm2 + 2*Vm3) - this->RR()*T*T05*Vm2*bpVm2)
        /(T05*Vm2*bpVm2*bmVm2);
}


//Real deviative dp/dT at constant molar volume
//(molar values)
inline scalar redlichKwong::dpdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return 0.5*a_/(pow(T, 1.5)*Vm*(b_ + Vm)) - this->RR()/(b_ - Vm);
}


//Real  deviative dv/dT at constant pressure
//using implicit differentiation
//(molar values)
inline scalar redlichKwong::dvdT(const scalar rho, const scalar T) const
{
    return -this->dpdT(rho, T)/this->dpdv(rho, T);
}


//Real  deviative dv/dp at constant temperature
//(molar values)
inline scalar redlichKwong::dvdp(const scalar rho, const scalar T)  const
{
    return 1/this->dpdv(rho, T);
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

    return this->RR()*T*log(Vm - b_) + a_/(b_*sqrt(T))*log((b_ + Vm)/Vm);
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

    return this->RR()*log(Vm - b_) - a_/(2*b_*pow(T, 1.5))*log((b_ + Vm)/Vm);
}


//(molar values)
inline scalar redlichKwong::d2pdT2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;

    return -0.75*a_/(pow(T, 2.5)*Vm*(b_ + Vm));
}


//(molar values)
inline scalar redlichKwong::d2pdv2(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar Vm3 = Vm2*Vm;
    scalar T05 = pow(T, 0.5);
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bmVm3 = bmVm2*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;
    scalar bpVm3 = bpVm2*bpVm;

    return 2*(a_*bmVm3*(b2_ + 3*b_*Vm + 3*Vm2) + this->RR()*T*T05*Vm3*bpVm3)
        /(T05*Vm3*bpVm3*bmVm3);
}


//(molar values)
//using second order implicit differentiation
inline scalar redlichKwong::d2vdT2(const scalar rho, const scalar T) const
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
inline scalar redlichKwong::d2pdvdT(const scalar rho, const scalar T) const
{
    scalar Vm = this->W()/rho;
    scalar Vm2 = Vm*Vm;
    scalar T15_ = pow(T, 1.5);
    scalar bmVm = b_ - Vm;
    scalar bmVm2 = bmVm*bmVm;
    scalar bpVm = b_ + Vm;
    scalar bpVm2 = bpVm*bpVm;

    return -(0.5*a_*bmVm2*(b_ + 2*Vm) + this->RR()*T15_*Vm2*bpVm2)
        /(T15_*Vm2*bpVm2*bmVm2);
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
    scalar T25 = pow(T, 2.5);

    return 0.75*a_/(T25*b_)*log((b_ +  Vm)/Vm);
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
               "inline scalar redlichKwong::rho"
               "(const scalar p, const scalar T, const scalar rho0) const "
           )   << "Maximum number of iterations exceeded"
               << abort(FatalError);
       }
   }
   while(mag(molarVolumePrevIteration-molarVolume) > tol_*(this->W()/rho0));

   return this->W()/molarVolume;
}


//- Return density [kg/m^3]
inline scalar redlichKwong::rho(const scalar p, const scalar T) const
{
    // using perfect gas equation as starting point
    return rho(p, T, p/(this->R()*T));
}


//- Return compressibility drho/dp at T=constant [s^2/m^2]
inline scalar redlichKwong::psi(const scalar rho, const scalar T) const
{
    return -this->dvdp(rho, T)*rho*rho/this->W();
}


//- Return compression factor []
inline scalar redlichKwong::Z
(
    const scalar p,
    const scalar T,
    const scalar rho0
) const
{
    return p/(this->R()*T*this->rho(p, T, rho0));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void redlichKwong::operator+=(const redlichKwong& rk)
{
    specie::operator+=(rk);
}


inline void redlichKwong::operator*=(const scalar s)
{
    specie::operator*=(s);
}


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
