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
    Foam::nasaHeatCapacityPolynomial

Description
    Nasa Heat Capacity Polynomial for perfect Gas (7.order polynominal) --> freely available for many fluids
    --> paper title: NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species

    templated into the equationOfState
    -> uses the equation of state to calculate all real Gas properties like Enthalpy, Entropy ...
    -> can not be used with the perfectGas equation of state

    Equations for the real gas correction: Have a look at thermodnamics books e.g. Thermodynamics: An Engineering Approch, 5 Edition, Chapter 12

SourceFiles
    nasaHeatCapacityPolynomialI.H
    nasaHeatCapacityPolynomial.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#ifndef nasaHeatCapacityPolynomial_H
#define nasaHeatCapacityPolynomial_H



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class nasaHeatCapacityPolynomial;

template<class equationOfState>
inline nasaHeatCapacityPolynomial<equationOfState> operator+
(
    const nasaHeatCapacityPolynomial<equationOfState>&,
    const nasaHeatCapacityPolynomial<equationOfState>&
);

template<class equationOfState>
inline nasaHeatCapacityPolynomial<equationOfState> operator-
(
    const nasaHeatCapacityPolynomial<equationOfState>&,
    const nasaHeatCapacityPolynomial<equationOfState>&
);

template<class equationOfState>
inline nasaHeatCapacityPolynomial<equationOfState> operator*
(
    const scalar,
    const nasaHeatCapacityPolynomial<equationOfState>&
);

template<class equationOfState>
inline nasaHeatCapacityPolynomial<equationOfState> operator==
(
    const nasaHeatCapacityPolynomial<equationOfState>&,
    const nasaHeatCapacityPolynomial<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const nasaHeatCapacityPolynomial<equationOfState>&
);


/*---------------------------------------------------------------------------*\
           Class nasaHeatCapacityPolynomial Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class nasaHeatCapacityPolynomial
:
    public equationOfState
{
    // Private data

        scalar a1_;
        scalar a2_;
        scalar a3_;
        scalar a4_;
        scalar a5_;
        scalar a6_;
        scalar a7_;
        scalar e0_std;
        scalar s0_std;
        scalar integral_p_dv_std;
        scalar integral_dpdT_dv_std;
        scalar cp_std;

    // Private member functions

        //- Construct from components
        //CL: used for the operator+
        inline nasaHeatCapacityPolynomial
        (
            const equationOfState& st,
            const scalar a1,
            const scalar a2,
            const scalar a3,
            const scalar a4,
            const scalar a5,
            const scalar a6,
            const scalar a7
    );

        //- Construct from components
        //CL: used for the operator*
        inline nasaHeatCapacityPolynomial
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
        );

public:

//Variable

    // Constructors
        //- Construct from Istream
        nasaHeatCapacityPolynomial(Istream&);

        //- Construct from dictionary
        nasaHeatCapacityPolynomial(const dictionary& dict);

        //- Construct as named copy
        inline nasaHeatCapacityPolynomial(const word&, const nasaHeatCapacityPolynomial&);

        //- Construct and return a clone
        inline autoPtr<nasaHeatCapacityPolynomial> clone() const;

        //- Selector from Istream
        inline static autoPtr<nasaHeatCapacityPolynomial> New(Istream& is);

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

        inline void operator+=(const nasaHeatCapacityPolynomial&);
        inline void operator-=(const nasaHeatCapacityPolynomial&);


    // Friend operators

        friend nasaHeatCapacityPolynomial operator+ <equationOfState>
        (
            const nasaHeatCapacityPolynomial&,
            const nasaHeatCapacityPolynomial&
        );

        friend nasaHeatCapacityPolynomial operator- <equationOfState>
        (
            const nasaHeatCapacityPolynomial&,
            const nasaHeatCapacityPolynomial&
        );

        friend nasaHeatCapacityPolynomial operator* <equationOfState>
        (
            const scalar,
            const nasaHeatCapacityPolynomial&
        );

        friend nasaHeatCapacityPolynomial operator== <equationOfState>
        (
            const nasaHeatCapacityPolynomial&,
            const nasaHeatCapacityPolynomial&
        );


    // IOstream Operators

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const nasaHeatCapacityPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nasaHeatCapacityPolynomialI.H"

#ifdef NoRepository
#   include "nasaHeatCapacityPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
