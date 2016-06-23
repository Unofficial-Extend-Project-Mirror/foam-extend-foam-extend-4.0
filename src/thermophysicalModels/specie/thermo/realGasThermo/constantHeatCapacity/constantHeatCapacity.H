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

Class
    Foam::constantHeatCapacity

Description
    real gas thermodynamic class --> constant perfect gas heat capacity
    Important: the perfect gas heat capacity is constant, the real heat capacity is not constant due to real gas corrections

    templated into the equationOfState
    -> uses the equation of state to calculate all real Gas properties like Enthalpy, Entropy ...
    -> can not be used with the perfectGas equation of state

    Equations for the real gas correction: Have a look at thermodnamics books e.g. Thermodynamics:
    An Engineering Approch, 5 Edition, Chapter 12

SourceFiles
    constantHeatCapacityI.H
    constantHeatCapacity.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#ifndef constantHeatCapacity_H
#define constantHeatCapacity_H



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class constantHeatCapacity;

template<class equationOfState>
inline constantHeatCapacity<equationOfState> operator+
(
    const constantHeatCapacity<equationOfState>&,
    const constantHeatCapacity<equationOfState>&
);

template<class equationOfState>
inline constantHeatCapacity<equationOfState> operator-
(
    const constantHeatCapacity<equationOfState>&,
    const constantHeatCapacity<equationOfState>&
);

template<class equationOfState>
inline constantHeatCapacity<equationOfState> operator*
(
    const scalar,
    const constantHeatCapacity<equationOfState>&
);

template<class equationOfState>
inline constantHeatCapacity<equationOfState> operator==
(
    const constantHeatCapacity<equationOfState>&,
    const constantHeatCapacity<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const constantHeatCapacity<equationOfState>&
);


/*---------------------------------------------------------------------------*\
           Class constantHeatCapacity Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class constantHeatCapacity
:
    public equationOfState
{
    // Private data

        //CL: specific values
        scalar Cp0_;
        //CL: molar values
        scalar cp0_;
        scalar e0_std;
        scalar s0_std;
        scalar integral_p_dv_std;
        scalar integral_dpdT_dv_std;
        scalar cp_std;

    // Private member functions

        //- Construct from components
        //CL: used for the operator+
        inline constantHeatCapacity
        (
            const equationOfState& st,
            const scalar cp0_
        );

        //- Construct from components
        //CL: used for the operator*
        inline constantHeatCapacity
        (
            const equationOfState& st,
            const scalar cp0_,
            const scalar e0_std_,
            const scalar s0_std_,
            const scalar integral_p_dv_std_,
            const scalar integral_dpdT_dv_std_,
            const scalar cp_std_
        );

public:

//Variable

    // Constructors
        //- Construct from Istream
        constantHeatCapacity(Istream&);

        //- Construct as named copy
        inline constantHeatCapacity(const word&, const constantHeatCapacity&);

        //- Construct and return a clone
        inline autoPtr<constantHeatCapacity> clone() const;

        //- Selector from Istream
        inline static autoPtr<constantHeatCapacity> New(Istream& is);

    // Member Functions

    //- perfect Gas Enthalpy [J/kmol]
        inline scalar h0(const scalar T) const;

        //- perfect Gas Entropy [J/(kmol K)]
        inline scalar s0(const scalar T) const;

        //- perfect Gas internal Energy  [J/kmol]
        inline scalar e0(const scalar T) const;

        //- perfect gas Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cv0(const scalar T) const;

        //- perfect gas Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp0(const scalar T) const;

        //- Limited Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp(const scalar rho, const scalar T) const;

        //- non Limited Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp_nonLimited(const scalar rho, const scalar T) const;

        //- Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cv(const scalar rho, const scalar T) const;

        //- Enthalpy [J/kmol]
        inline scalar h(const scalar rho, const scalar T) const;

        //- Entropy [J/(kmol K)]
        inline scalar s(const scalar rho,const scalar T) const;

        //- Internal Energy [J/kmol]
        inline scalar e(const scalar rho, const scalar T) const;


    // Member operators

        inline void operator+=(const constantHeatCapacity&);
        inline void operator-=(const constantHeatCapacity&);


    // Friend operators

        friend constantHeatCapacity operator+ <equationOfState>
        (
            const constantHeatCapacity&,
            const constantHeatCapacity&
        );

        friend constantHeatCapacity operator- <equationOfState>
        (
            const constantHeatCapacity&,
            const constantHeatCapacity&
        );

        friend constantHeatCapacity operator* <equationOfState>
        (
            const scalar,
            const constantHeatCapacity&
        );

        friend constantHeatCapacity operator== <equationOfState>
        (
            const constantHeatCapacity&,
            const constantHeatCapacity&
        );


    // IOstream Operators

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const constantHeatCapacity&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constantHeatCapacityI.H"

#ifdef NoRepository
#   include "constantHeatCapacity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
