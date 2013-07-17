/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Class
    Foam::viewFactorRadiation

Description
    helper class to sum up external radiation sources

      
SourceFiles
    viewFactorRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef viewFactorRadiation_H
#define viewFactorRadiation_H

#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class viewFactorRadiation;
Ostream& operator<<(Ostream&, const viewFactorRadiation&);

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class viewFactorRadiation
{
    // Private data

        //- name
        word name_;

        //- Tinf [K]
        scalar Tinf_;

        //- view factor [-]
        scalarField F_;

        //- emissivity [-]
        scalar epsilon_;


public:

    // Selectors

        //- Return a reference to the selected thermal model
        static autoPtr<viewFactorRadiation> New
        (
            const word& name,
            const dictionary& dict,
            const label size
        )
        {
            return autoPtr<viewFactorRadiation>
            (
                new viewFactorRadiation(name, dict, size)
            );
        };

    // Constructors

        //- Construct from components
        explicit viewFactorRadiation
        (
            const word& name,
            const dictionary& dict,
            const label size
        );


        //- Construct and return a clone
        Foam::autoPtr<viewFactorRadiation> clone() const
        {
            return autoPtr<viewFactorRadiation>
            (
                new viewFactorRadiation(*this)
            );
        }


    // Member functions

        // Access

        // Mapping functions

        // Evaluation functions

            //- Return Tinf
            scalar Tinf() const
            {
                return Tinf_;
            };

            //- Return F
            const scalarField& F() const
            {
                return F_;
            };

            //- Return epsilon
            scalar epsilon() const
            {
                return epsilon_;
            };


        //- Write
        void write(Ostream&) const;

    // IOstream operators

        friend Ostream& operator<<(Ostream&, const viewFactorRadiation&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
