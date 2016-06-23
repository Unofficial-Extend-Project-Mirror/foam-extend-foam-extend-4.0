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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
//CL: used for the operator+
template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState>::constantHeatCapacity
(
    const equationOfState& st,
    const scalar cp0_
)
:
    equationOfState(st),
    cp0_(cp0_),
    e0_std(e0(this->Tstd)),
    s0_std(s0(this->Tstd)),
    integral_p_dv_std(this->integral_p_dv(this->rhostd(),this->Tstd)),
    integral_dpdT_dv_std(this->integral_dpdT_dv(this->rhostd(),this->Tstd)),
    cp_std(this->cp_nonLimited(this->rhostd(),this->Tstd))
{}


//- Construct from components
//CL: used for the operator*
template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState>::constantHeatCapacity
(
    const equationOfState& st,
    const scalar cp0_,
    const scalar e0_std_,
    const scalar s0_std_,
    const scalar integral_p_dv_std_,
    const scalar integral_dpdT_dv_std_,
    const scalar cp_std_
)
:
    equationOfState(st),
    cp0_(cp0_),
    e0_std(e0_std_),
    s0_std(s0_std_),
    integral_p_dv_std(integral_p_dv_std_),
    integral_dpdT_dv_std(integral_dpdT_dv_std_),
    cp_std(cp_std_)
{}


template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState>::constantHeatCapacity
(
    const word& name,
    const constantHeatCapacity& np
)
:
    equationOfState(name, np),
    cp0_(np.cp0_),
    e0_std(np.e0_std),
    s0_std(np.s0_std),
    integral_p_dv_std(np.integral_p_dv_std),
    integral_dpdT_dv_std(np.integral_dpdT_dv_std),
    cp_std(np.cp_std)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::constantHeatCapacity<equationOfState> >
Foam::constantHeatCapacity<equationOfState>::clone() const
{
    return autoPtr<constantHeatCapacity<equationOfState> >
    (
        new constantHeatCapacity<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::constantHeatCapacity<equationOfState> >
Foam::constantHeatCapacity<equationOfState>::New(Istream& is)
{
    return autoPtr<constantHeatCapacity<equationOfState> >
    (
        new constantHeatCapacity<equationOfState>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//used to calculate the internal energy
//perfect gas enthalpy
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::h0
(
    const scalar T
) const
{
    return cp0_*T;
}


//used to calculate the internal energy
//perfect gas internal energy
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::e0
(
    const scalar T
) const
{
    return this->h0(T) - this->RR()*T;
}


// used to calculate the entropy
// perfect gas entropy
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::s0
(
    const scalar T
) const
{
     return cp0_*log(T);
}


//perfect gas cp
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::cp0
(
    const scalar T
) const
{
    return cp0_;
}


//perfect gas cv
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::cv0
(
    const scalar T
) const
{
    return this->cp0(T)-this->RR();
}


//function to calculate real gas cp
//using cp=cv+(dp/dT)^2/(dp/dv)
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::cp
(
     const scalar rho,
     const scalar T
) const
{
    // CL: within the vapor dome, some stabilty problems might occur due to unphysical results (EOS not valite in vapor dome)
    // CL: In order to find the right solution, the newton solver might go into the vapor dome (in an intermediate step)
    // CL: Therefore, to avoid a crash, the following stabilisation step should be used
    // CL: dpdv(rho,T)  is =0 at some points within the vapour dome. To increase stability, (dp/dv) has to be limited
    // CL: cp can be negative within the vapor dome. To avoid this nonphysical result, the absolute value is used.
    // CL: within the vapourdome and at the critical point, cp increases to very high values --> infinity,
    // CL: this would decrease the stability, so cp will be limited to 20 times the cp @ STD

    return
    min
    (
        cp_std*20,
        fabs
        (
            this->cv(rho,T)
          - T*pow((this->dpdT(rho, T)), 2)/min(this->dpdv(rho, T), -1)
        )
    );
}


// this function is needed to get cp @ STD (without the limit imposed in the function above),
// which in turn is needed to limit the cp in the function above
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::cp_nonLimited
(
     const scalar rho,
     const scalar T
) const
{
    return
        fabs(this->cv(rho, T)
      - T*pow((this->dpdT(rho, T)), 2)/min(this->dpdv(rho, T), -1));
}


//function to calculate real gas c
//cv=cv0+T*integral d2p/dT2 dv
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::cv
(
    const scalar rho,
    const scalar T
) const
{
    return this->cv0(T) + T*this->integral_d2pdT2_dv(rho, T);
}


//function to calculate real gas enthalpy
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::h
(
    const scalar rho,
    const scalar T
) const
{
    return
        this->e(rho, T)
      + this->p(rho, T)/rho*this->W()
      - this->Pstd()/this->rhostd()*this->W();
}


// function to calculate real gas internal energy
// important assumption used: internal Energie is 0 at STD conditions.
// equation: du= cv0 dT +[T*dp/dT -p]dv
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::e
(
    const scalar rho,
    const scalar T
) const
{
    return
      - this->Tstd()*integral_dpdT_dv_std
      + integral_p_dv_std
      + this->e0(T)
      - e0_std
      + T*this->integral_dpdT_dv(rho, T)
      - this->integral_p_dv(rho, T);
}


//function to calculate real gas entropy
// important assumption used: the Entropy is 0 at STD conditions.
// equation: ds= cv0/T * dT +  dp/dT *dv
// --> integral cv0/T dT = s0(T) -s0(Tstd) - R*ln(T/Tstd)  --> due to s0(T)-s0(Tstd)=integral cp0/T dT
template<class equationOfState>
inline Foam::scalar Foam::constantHeatCapacity<equationOfState>::s
(
    const scalar rho,
    const scalar T
) const
{
    return
      - integral_dpdT_dv_std
      + this->s0(T)
      - s0_std
      - this->RR()*log(T/this->Tstd())
      + this->integral_dpdT_dv(rho, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::constantHeatCapacity<equationOfState>::operator+=
(
    const constantHeatCapacity<equationOfState>& np
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(np);

    molr1 /= this->nMoles();
    scalar molr2 = np.nMoles()/this->nMoles();

    cp0_ = molr1*cp0_ + molr2*np.cp0_;
}


template<class equationOfState>
inline void Foam::constantHeatCapacity<equationOfState>::operator-=
(
    const constantHeatCapacity<equationOfState>& np
)
{
    scalar molr1 = this->nMoles();

    constantHeatCapacity::operator-=(np);

    molr1 /= this->nMoles();
    scalar molr2 = np.nMoles()/this->nMoles();

    cp0_ = molr1*cp0_ - molr2*np.cp0_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState> Foam::operator+
(
    const constantHeatCapacity<equationOfState>& np1,
    const constantHeatCapacity<equationOfState>& np2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(np1)
      + static_cast<const equationOfState&>(np2)
    );

    //CL: Important, calls a different constructor as operator*
    //CL: the coefficients as well as the EOS (coefficients) changed
    //CL: therefore, the values at STD needs to be recalculated
    return constantHeatCapacity<equationOfState>
    (
        eofs,
        np1.nMoles()/eofs.nMoles()*np1.cp0_
      + np2.nMoles()/eofs.nMoles()*np2.cp0_
    );
}


template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState> Foam::operator-
(
    const constantHeatCapacity<equationOfState>& np1,
    const constantHeatCapacity<equationOfState>& np2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(np1)
      - static_cast<const equationOfState&>(np2)
    );

    return constantHeatCapacity<equationOfState>
    (
        eofs,
        np1.nMoles()/eofs.nMoles()*np1.cp0_
      - np2.nMoles()/eofs.nMoles()*np2.cp0_
    );
}


template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState> Foam::operator*
(
    const scalar s,
    const constantHeatCapacity<equationOfState>& np
)
{
    //CL: values at STD don't need to be recalculated,
    //CL: therefore, providing the values in the constructor
    return constantHeatCapacity<equationOfState>
    (
        s*static_cast<const equationOfState&>(np),
        np.cp0_,
        np.e0_std,
        np.s0_std,
        np.integral_p_dv_std,
        np.integral_dpdT_dv_std,
        np.cp_std
    );
}


template<class equationOfState>
inline Foam::constantHeatCapacity<equationOfState> Foam::operator==
(
    const constantHeatCapacity<equationOfState>& np1,
    const constantHeatCapacity<equationOfState>& np2
)
{
    return np2 - np1;
}


// ************************************************************************* //
