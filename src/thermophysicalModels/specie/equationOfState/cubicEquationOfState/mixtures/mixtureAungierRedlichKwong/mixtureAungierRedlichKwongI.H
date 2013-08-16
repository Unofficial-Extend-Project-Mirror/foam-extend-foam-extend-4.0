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

#include "mixtureAungierRedlichKwong.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline mixtureAungierRedlichKwong::mixtureAungierRedlichKwong
(
    const aungierRedlichKwong& sp,
    scalar Tcrit,
    scalar azentricFactor,
    scalar Vcrit,
    scalar Zcrit,
    scalar rhoMin,
    scalar rhoMax
)
:
    aungierRedlichKwong(sp),    
    Vcrit_(Vcrit),
    Zcrit_(Zcrit)
{
    //CL: Saving critical data
    Tcrit_=Tcrit;
    pcrit_=Zcrit*this->RR*Tcrit/Vcrit;  
    rhocrit_=this->W()/Vcrit_;    
    azentricFactor_=azentricFactor;

    //CL: calculating the aungier redlich kwong coefficience
    a0_=0.42747*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_;
    b_=0.08664*this->RR*Tcrit_/pcrit_;
    c_=this->RR*Tcrit_/(pcrit_+(a0_/(this->W()/rhocrit_*(this->W()/rhocrit_+b_))))+b_-this->W()/rhocrit_;
    n_=0.4986+1.2735*azentricFactor_+0.4754*pow(azentricFactor_,2);

    //CL: set rhoMax and rhoMin
    rhoMin_=rhoMin;
    rhoMax_=rhoMax;

    //CL: 
    b2_=pow(b_,2);
    b3_=pow(b_,3);
    b4_=pow(b_,4);
    b5_=pow(b_,5);
    c2_=pow(c_,2);

    //CL: Starting GUESS for the density by ideal gas law
    rhostd_=this->rho(Pstd,Tstd,Pstd/(Tstd*this->R()));
}




// Construct as named copy
inline mixtureAungierRedlichKwong::mixtureAungierRedlichKwong(const word& name, const mixtureAungierRedlichKwong& pg)
:
    aungierRedlichKwong(name, pg),
    Vcrit_(pg.Vcrit_),
    Zcrit_(pg.Zcrit_)
{}



// Construct and return a clone
inline autoPtr<mixtureAungierRedlichKwong> mixtureAungierRedlichKwong::clone() const
{
    return autoPtr<mixtureAungierRedlichKwong>(new mixtureAungierRedlichKwong(*this));
}


// Selector from Istream
inline autoPtr<mixtureAungierRedlichKwong> mixtureAungierRedlichKwong::New(Istream& is)
{
    return autoPtr<mixtureAungierRedlichKwong>(new mixtureAungierRedlichKwong(is));
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline void mixtureAungierRedlichKwong::operator+=(const mixtureAungierRedlichKwong& pr)
{
    scalar molr1 = this->nMoles();

    aungierRedlichKwong::operator+=(pr);

    molr1 /= this->nMoles();
    scalar molr2 = pr.nMoles()/this->nMoles();

    //CL: calculating new critical point data
    Tcrit_=molr1*Tcrit_+molr2*pr.Tcrit_;
    Zcrit_=molr1*Zcrit_+molr2*pr.Zcrit_;
    Vcrit_=molr1*Vcrit_+molr2*pr.Vcrit_;
    pcrit_=Zcrit_*this->RR*Tcrit_/Vcrit_;
    rhocrit_=this->W()/Vcrit_;
    
    //CL: calculating new azentric factor
    azentricFactor_=molr1*azentricFactor_+molr2*pr.azentricFactor_;

    //CL: set rhoMax and rhoMin
    rhoMin_=min(rhoMin_, pr.rhoMin_);
    rhoMax_=max(rhoMax_, pr.rhoMax_);

    //CL: calculating the aungier redlich kwong coefficience
    a0_=0.42747*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_;
    b_=0.08664*this->RR*Tcrit_/pcrit_;
    c_=this->RR*Tcrit_/(pcrit_+(a0_/(this->W()/rhocrit_*(this->W()/rhocrit_+b_))))+b_-this->W()/rhocrit_;
    n_=0.4986+1.2735*azentricFactor_+0.4754*pow(azentricFactor_,2);

    //CL: 
    b2_=pow(b_,2);
    b3_=pow(b_,3);
    b4_=pow(b_,4);
    b5_=pow(b_,5);
    c2_=pow(c_,2);

    //CL: Starting GUESS for the density by ideal gas law
    rhostd_=this->rho(Pstd,Tstd,Pstd/(Tstd*this->R()));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline mixtureAungierRedlichKwong operator+
(
    const mixtureAungierRedlichKwong& pr1,
    const mixtureAungierRedlichKwong& pr2
)
{

    aungierRedlichKwong ark
    (
        static_cast<const aungierRedlichKwong&>(pr1)
       +static_cast<const aungierRedlichKwong&>(pr2)
    );

    return mixtureAungierRedlichKwong
           (
               ark,
               pr1.nMoles()/ark.nMoles()*pr1.Tcrit_
              +pr2.nMoles()/ark.nMoles()*pr2.Tcrit_,
               pr1.nMoles()/ark.nMoles()*pr1.azentricFactor_
              +pr2.nMoles()/ark.nMoles()*pr2.azentricFactor_,
               pr1.nMoles()/ark.nMoles()*pr1.Vcrit_
              +pr2.nMoles()/ark.nMoles()*pr2.Vcrit_,
               pr1.nMoles()/ark.nMoles()*pr1.Zcrit_
              +pr2.nMoles()/ark.nMoles()*pr2.Zcrit_,
               min(pr1.rhoMin_,pr2.rhoMin_),
               max(pr1.rhoMax_,pr2.rhoMax_)
           );
}


inline mixtureAungierRedlichKwong operator*
(
    const scalar s,
    const mixtureAungierRedlichKwong& pr
)
{
    return mixtureAungierRedlichKwong
           (
               s*static_cast<const aungierRedlichKwong&>(pr),
               pr.Tcrit_, 
               pr.azentricFactor_, 
               pr.Vcrit_, 
               pr.Zcrit_, 
               pr.rhoMin_, 
               pr.rhoMax_
           );
}




} // End namespace Foam

// ************************************************************************* //
