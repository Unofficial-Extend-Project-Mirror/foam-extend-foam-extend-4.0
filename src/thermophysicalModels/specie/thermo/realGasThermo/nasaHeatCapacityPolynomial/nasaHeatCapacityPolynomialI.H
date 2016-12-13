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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
//CL: used for the operator+
template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState>::nasaHeatCapacityPolynomial
(
    const equationOfState& st,
    const scalar a1,
    const scalar a2,
    const scalar a3,
    const scalar a4,
    const scalar a5,
    const scalar a6,
    const scalar a7
)
:
    equationOfState(st),
    a1_(a1),
    a2_(a2),
    a3_(a3),
    a4_(a4),
    a5_(a5),
    a6_(a6),
    a7_(a7),
    e0_std(e0(this->Tstd())),
    s0_std(s0(this->Tstd())),
    integral_p_dv_std(this->integral_p_dv(this->rhostd(),this->Tstd())),
    integral_dpdT_dv_std(this->integral_dpdT_dv(this->rhostd(),this->Tstd())),
    cp_std(this->cp_nonLimited(this->rhostd(),this->Tstd()))
{}

//- Construct from components
//CL: used for the operator*
template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState>::nasaHeatCapacityPolynomial
(
    const equationOfState& st,
    const scalar a1,
    const scalar a2,
    const scalar a3,
    const scalar a4,
    const scalar a5,
    const scalar a6,
    const scalar a7,
    const scalar e0_std_,
    const scalar s0_std_,
    const scalar integral_p_dv_std_,
    const scalar integral_dpdT_dv_std_,
    const scalar cp_std_
)
:
    equationOfState(st),
    a1_(a1),
    a2_(a2),
    a3_(a3),
    a4_(a4),
    a5_(a5),
    a6_(a6),
    a7_(a7),
    e0_std(e0_std_),
    s0_std(s0_std_),
    integral_p_dv_std(integral_p_dv_std_),
    integral_dpdT_dv_std(integral_dpdT_dv_std_),
    cp_std(cp_std_)
{}



template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState>::nasaHeatCapacityPolynomial
(
    const word& name,
    const nasaHeatCapacityPolynomial& np
)
:
    equationOfState(name, np),
    a1_(np.a1_),
    a2_(np.a2_),
    a3_(np.a3_),
    a4_(np.a4_),
    a5_(np.a5_),
    a6_(np.a6_),
    a7_(np.a7_),
    e0_std(np.e0_std),
    s0_std(np.s0_std),
    integral_p_dv_std(np.integral_p_dv_std),
    integral_dpdT_dv_std(np.integral_dpdT_dv_std),
    cp_std(np.cp_std)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::nasaHeatCapacityPolynomial<equationOfState> >
Foam::nasaHeatCapacityPolynomial<equationOfState>::clone() const
{
    return autoPtr<nasaHeatCapacityPolynomial<equationOfState> >
    (
        new nasaHeatCapacityPolynomial<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::nasaHeatCapacityPolynomial<equationOfState> >
Foam::nasaHeatCapacityPolynomial<equationOfState>::New(Istream& is)
{
    return autoPtr<nasaHeatCapacityPolynomial<equationOfState> >
    (
        new nasaHeatCapacityPolynomial<equationOfState>(is)
    );
}


 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//used to calculate the internal energy
//perfect gas enthalpy
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::h0
(
    const scalar T
) const
{
    return
    this->RR()*
    (
        - this->a1_/T
        + this->a2_*log(T)
        + T*
          (
              this->a3_ + T*
              (
                  0.5*this->a4_ + T*
                  (
                      this->a5_/3 + T*
                      (
                          0.25*this->a6_ + T*(0.2*this->a7_)
                      )
                  )
              )
          )
    );
}


//used to calculate the internal energy
//perfect gas internal energy
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::e0
(
    const scalar T
) const
{
    return this->h0(T) - this->RR()*T;
}
// used to calculate the entropy
// perfect gas entropy
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::s0
(
    const scalar T
) const
{
    scalar oneOverT = 1/T;

    return this->RR()*
    (
        oneOverT*(- 0.5*this->a1_*oneOverT - this->a2_)
      + this->a3_*log(T)
      + T*
        (
            this->a4_ + T*
            (
                0.5*this->a5_ + T*
                (
                    this->a6_/3 + T*(0.25*this->a7_)
                )
            )
        )
    );
}


//perfect gas cp
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::cp0
(
    const scalar T
) const
{
    scalar oneOverT = 1/T;

    return this->RR()*
    (
        oneOverT*(this->a1_*oneOverT + this->a2_)
      + this->a3_
      + T*(this->a4_ + T*(this->a5_ + T*(this->a6_ + T*this->a7_)))
    );
}


//perfect gas cv
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::cv0
(
    const scalar T
) const
{

    return this->cp0(T)-this->RR();
}


//function to calculate real gas cp
//using cp=cv+(dp/dT)^2/(dp/dv)
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::cp
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

    return min
    (
        cp_std*20,
        fabs
        (
            this->cv(rho, T)
          - T*pow(this->dpdT(rho, T), 2)/min(this->dpdv(rho, T), -1)
        )
    );
}

// this function is needed to get cp @ STD (without the limit imposed in the function above),
// which in turn is needed to limit the cp in the function above
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::cp_nonLimited
(
     const scalar rho,
     const scalar T
) const
{
     return fabs
     (
         this->cv(rho,T)-T*pow((this->dpdT(rho, T)),2)/min(this->dpdv(rho, T), -1)
     );
}


//function to calculate real gas cv
//cv=cv0+T*integral d2p/dT2 dv
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::cv
(
    const scalar rho,
    const scalar T
) const
{
    return this->cv0(T)+T*this->integral_d2pdT2_dv(rho, T);
}


//function to calculate real gas enthalpy
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::h
(
    const scalar rho,
    const scalar T
) const
{
    return this->e(rho,T)
      + this->p(rho, T)/rho*this->W()
      - this->Pstd()/this->rhostd()*this->W();
}


// function to calculate real gas internal energy
// important assumption used: internal Energie is 0 at STD conditions.
// equation: du= cv0 dT +[T*dp/dT -p]dv
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::e
(
    const scalar rho,
    const scalar T
) const
{
    return
    (
      - this->Tstd()*integral_dpdT_dv_std
      + integral_p_dv_std
      + this->e0(T)
      - e0_std
      + T*this->integral_dpdT_dv(rho, T)
      - this->integral_p_dv(rho, T)
    );
}


//function to calculate real gas entropy
// important assumption used: the Entropy is 0 at STD conditions.
// equation: ds= cv0/T * dT +  dp/dT *dv
// --> integral cv0/T dT = s0(T) -s0(Tstd) - R*ln(T/Tstd)  --> due to s0(T)-s0(Tstd)=integral cp0/T dT
template<class equationOfState>
inline Foam::scalar Foam::nasaHeatCapacityPolynomial<equationOfState>::s
(
    const scalar rho,
    const scalar T
) const
{
    return
    (
      - integral_dpdT_dv_std
      + this->s0(T)
      - s0_std
      - this->RR()*log(T/this->Tstd())
      + this->integral_dpdT_dv(rho, T)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::nasaHeatCapacityPolynomial<equationOfState>::operator+=
(
    const nasaHeatCapacityPolynomial<equationOfState>& np
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(np);

    molr1 /= this->nMoles();
    scalar molr2 = np.nMoles()/this->nMoles();

    a1_ = molr1*a1_ + molr2*np.a1_;
    a2_ = molr1*a2_ + molr2*np.a2_;
    a3_ = molr1*a3_ + molr2*np.a3_;
    a4_ = molr1*a4_ + molr2*np.a4_;
    a5_ = molr1*a5_ + molr2*np.a5_;
    a6_ = molr1*a6_ + molr2*np.a6_;
    a7_ = molr1*a7_ + molr2*np.a7_;
}


template<class equationOfState>
inline void Foam::nasaHeatCapacityPolynomial<equationOfState>::operator-=
(
    const nasaHeatCapacityPolynomial<equationOfState>& np
)
{
    scalar molr1 = this->nMoles();

    nasaHeatCapacityPolynomial::operator-=(np);

    molr1 /= this->nMoles();
    scalar molr2 = np.nMoles()/this->nMoles();

    a1_ = molr1*a1_ - molr2*np.a1_;
    a2_ = molr1*a2_ - molr2*np.a2_;
    a3_ = molr1*a3_ - molr2*np.a3_;
    a4_ = molr1*a4_ - molr2*np.a4_;
    a5_ = molr1*a5_ - molr2*np.a5_;
    a6_ = molr1*a6_ - molr2*np.a6_;
    a7_ = molr1*a7_ - molr2*np.a7_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState> Foam::operator+
(
    const nasaHeatCapacityPolynomial<equationOfState>& np1,
    const nasaHeatCapacityPolynomial<equationOfState>& np2
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
    return nasaHeatCapacityPolynomial<equationOfState>
    (
        eofs,
        np1.nMoles()/eofs.nMoles()*np1.a1_
      + np2.nMoles()/eofs.nMoles()*np2.a1_,
        np1.nMoles()/eofs.nMoles()*np1.a2_
      + np2.nMoles()/eofs.nMoles()*np2.a2_,
        np1.nMoles()/eofs.nMoles()*np1.a3_
      + np2.nMoles()/eofs.nMoles()*np2.a3_,
        np1.nMoles()/eofs.nMoles()*np1.a4_
      + np2.nMoles()/eofs.nMoles()*np2.a4_,
        np1.nMoles()/eofs.nMoles()*np1.a5_
      + np2.nMoles()/eofs.nMoles()*np2.a5_,
        np1.nMoles()/eofs.nMoles()*np1.a6_
      + np2.nMoles()/eofs.nMoles()*np2.a6_,
        np1.nMoles()/eofs.nMoles()*np1.a7_
      + np2.nMoles()/eofs.nMoles()*np2.a7_
    );
}


template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState> Foam::operator-
(
    const nasaHeatCapacityPolynomial<equationOfState>& np1,
    const nasaHeatCapacityPolynomial<equationOfState>& np2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(np1)
      - static_cast<const equationOfState&>(np2)
    );

    return nasaHeatCapacityPolynomial<equationOfState>
    (
        eofs,
        np1.nMoles()/eofs.nMoles()*np1.a1_
      - np2.nMoles()/eofs.nMoles()*np2.a1_,
        np1.nMoles()/eofs.nMoles()*np1.a2_
      - np2.nMoles()/eofs.nMoles()*np2.a2_,
        np1.nMoles()/eofs.nMoles()*np1.a3_
      - np2.nMoles()/eofs.nMoles()*np2.a3_,
        np1.nMoles()/eofs.nMoles()*np1.a4_
      - np2.nMoles()/eofs.nMoles()*np2.a4_,
        np1.nMoles()/eofs.nMoles()*np1.a5_
      - np2.nMoles()/eofs.nMoles()*np2.a5_,
        np1.nMoles()/eofs.nMoles()*np1.a6_
      - np2.nMoles()/eofs.nMoles()*np2.a6_,
        np1.nMoles()/eofs.nMoles()*np1.a7_
      - np2.nMoles()/eofs.nMoles()*np2.a7_
    );
}


template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState> Foam::operator*
(
    const scalar s,
    const nasaHeatCapacityPolynomial<equationOfState>& np
)
{
    //CL: values at STD don't need to be recalculated,
    //CL: therefore, providing the values in the constructor
    return nasaHeatCapacityPolynomial<equationOfState>
    (
        s*static_cast<const equationOfState&>(np),
        np.a1_,
        np.a2_,
        np.a3_,
        np.a4_,
        np.a5_,
        np.a6_,
        np.a7_,
        np.e0_std,
        np.s0_std,
        np.integral_p_dv_std,
        np.integral_dpdT_dv_std,
        np.cp_std
    );
}


template<class equationOfState>
inline Foam::nasaHeatCapacityPolynomial<equationOfState> Foam::operator==
(
    const nasaHeatCapacityPolynomial<equationOfState>& np1,
    const nasaHeatCapacityPolynomial<equationOfState>& np2
)
{
    return np2 - np1;
}


// ************************************************************************* //
