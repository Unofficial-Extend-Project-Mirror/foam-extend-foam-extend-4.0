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
    Foam::realGasSpecieThermo

Description

    This is an subversion of the specieThermo function. While the "orginal" specieThermo function
    is only valid for perfect Gases, this function is valid for real gases.

    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    realGasSpecieThermoI.H
    realGasSpecieThermo.C

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#ifndef realGasSpecieThermo_H
#define realGasSpecieThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class thermo> class realGasSpecieThermo;

template<class thermo>
inline realGasSpecieThermo<thermo> operator+
(
    const realGasSpecieThermo<thermo>&,
    const realGasSpecieThermo<thermo>&
);

template<class thermo>
inline realGasSpecieThermo<thermo> operator-
(
    const realGasSpecieThermo<thermo>&,
    const realGasSpecieThermo<thermo>&
);

template<class thermo>
inline realGasSpecieThermo<thermo> operator*
(
    const scalar,
    const realGasSpecieThermo<thermo>&
);

template<class thermo>
inline realGasSpecieThermo<thermo> operator==
(
    const realGasSpecieThermo<thermo>&,
    const realGasSpecieThermo<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const realGasSpecieThermo<thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class realGasSpecieThermo Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class realGasSpecieThermo
:
    public thermo
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const debug::tolerancesSwitch tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const debug::optimisationSwitch maxIter_;


    // Private member functions

        // return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(rho,T) and dF(rho,T)/dT
        inline void T
        (
            scalar f,
            scalar &T0,
            scalar p,
            scalar &rho0,
            scalar (realGasSpecieThermo::*F)(const scalar,const scalar) const,
            scalar (realGasSpecieThermo::*dFdT)(const scalar,const scalar) const
        ) const;

public:

        // Constructors

        //- construct from components
        inline realGasSpecieThermo(const thermo& sp);

        //- Construct from Istream
        realGasSpecieThermo(Istream&);

        //- Construct as named copy
        inline realGasSpecieThermo(const word& name, const realGasSpecieThermo&);


    // Member Functions

        // Fundamaental properties
        // (These functions must be provided in derived types)

            // Sensible enthalpy [J/kmol]
            //virtual scalar hs(const scalar) const;

            // Chemical enthalpy [J/kmol]
            //virtual scalar hc(const scalar) const;

        // Calculate and return derived properties
        // (These functions need not provided in derived types)

                //CL: isentropic expansion factor "gamma" (heat capacity ratio for perfect gas)
                inline scalar gamma(const scalar T, const scalar rho) const;

                //- Sensible internal energy [J/kmol]
                // inline scalar es(const scalar p, const scalar rho) const;

                //- Gibbs free energy [J/kmol]
                inline scalar g(const scalar rho,  const scalar T ) const;

                //- Helmholtz free energy [J/kmol]
                inline scalar a(const scalar rho,  const scalar T ) const;


            // Mass specific properties

                //- Heat capacity at constant pressure [J/(kg K)]
                inline scalar Cp(const scalar rho, const scalar T) const;

                //- Heat capacity at constant volume [J/(kg K)]
                inline scalar Cv(const scalar rho, const scalar T) const;

                //- Enthalpy [J/kg]
                inline scalar H(const scalar rho, const scalar T) const;

                //- Sensible enthalpy [J/kg]
               // inline scalar Hs(const scalar T) const;

                //- Chemical enthalpy [J/kg]
               // inline scalar Hc() const;

                //- Entropy [J/(kg K)]
                inline scalar S(const scalar rho, const scalar T) const;

                //- Internal energy [J/kg]
                inline scalar E(const scalar rho, const scalar T) const;

                //- Gibbs free energy [J/kg]
                inline scalar G(const scalar rho,  const scalar T) const;

                //- Helmholtz free energy [J/kg]
                inline scalar A(const scalar rho, const scalar T) const;

                //CL: Other variables

            //- Return compressibility drho/dp at h=constant [s^2/m^2]
                inline scalar psiH(const scalar rho, const scalar T) const;

                //- Return compressibility drho/dp at e=constant [s^2/m^2]
                inline scalar psiE(const scalar rho, const scalar T) const;

            //- Return compressibility drho/dH at p=constant
                inline scalar drhodH(const scalar rho, const scalar T) const;

            //- Return compressibility drho/dE at p=constant
                inline scalar drhodE(const scalar rho, const scalar T) const;

        // Energy->temperature  inversion functions

            //- Temperature from Enthalpy given an initial temperature T0
            inline void TH(const scalar H,  scalar &T0,const scalar p, scalar &psi0) const;

            //- Temperature from internal energy given an initial temperature T0
            inline void TE(const scalar E, scalar &T0,const scalar p, scalar &psi0) const;


    // Member operators

        inline void operator+=(const realGasSpecieThermo&);
        inline void operator-=(const realGasSpecieThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend realGasSpecieThermo operator+ <thermo>
        (
            const realGasSpecieThermo&,
            const realGasSpecieThermo&
        );

        friend realGasSpecieThermo operator- <thermo>
        (
            const realGasSpecieThermo&,
            const realGasSpecieThermo&
        );

        friend realGasSpecieThermo operator* <thermo>
        (
            const scalar s,
            const realGasSpecieThermo&
        );

        friend realGasSpecieThermo operator== <thermo>
        (
            const realGasSpecieThermo&,
            const realGasSpecieThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const realGasSpecieThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "realGasSpecieThermoI.H"

#ifdef NoRepository
#   include "realGasSpecieThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
