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
    Foam::viewFactor

Description
    helper class to sum up external radiation sources

      
SourceFiles
    viewFactor.C

\*---------------------------------------------------------------------------*/

#ifndef viewFactor_H
#define viewFactor_H

#include "externalRadiationSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class viewFactor
:
    public externalRadiationSource
{
    // Private data

        //- temperature at infinity
        scalar Tinf_;

        //- view factor
        scalarField F_;

        //- espison
        scalar epsilon_;

public:

    //- Runtime type information
    TypeName("viewFactor");


    // Constructors

        //- Construct from components
        explicit viewFactor
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );


        //- Construct from components
        explicit viewFactor
        (
            const word& name,
            const dictionary& dict
        );


        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            return autoPtr<externalRadiationSource>
            (
                new viewFactor(*this)
            );
        }

    // Destructor

        virtual ~viewFactor()
        {};


    // Member functions

        // Access

        const scalarField& F() const
        {
            return F_;
        };

        scalarField& F()
        {
            return F_;
        };

        const scalar& epsilon() const
        {
            return epsilon_;
        };

        scalar& epsilson()
        {
            return epsilon_;
        };


        // Mapping functions

        // Evaluation functions

            //- Return heat flux
            virtual tmp<scalarField> q(const scalarField& Tw) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
