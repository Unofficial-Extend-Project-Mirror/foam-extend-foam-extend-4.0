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
    Foam::pointSource

Description
    helper class to sum up external radiation sources

      
SourceFiles
    pointSource.C

\*---------------------------------------------------------------------------*/

#ifndef pointSource_H
#define pointSource_H

#include "constantFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pointSource
:
    public constantFlux
{
    // Private data

        //- Strenght of the point source
        scalar qmax_;

        //- radiation absorptivity
        scalar alpha_;

        //- Direction of the pointSource
        vector direction_;

public:

    //- Runtime type information
    TypeName("pointSource");


    // Constructors

        //- Construct from components
        explicit pointSource
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );


        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            return autoPtr<externalRadiationSource>
            (
                new pointSource(*this)
            );
        }

    // Destructor

        virtual ~pointSource()
        {};


    // Member functions

        // Access

        // Mapping functions

        // Evaluation functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
