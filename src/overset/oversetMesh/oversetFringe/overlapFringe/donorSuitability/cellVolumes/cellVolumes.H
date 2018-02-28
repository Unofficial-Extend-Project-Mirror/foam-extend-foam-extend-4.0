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

Class
    Foam::donorSuitability::cellVolumes

Description
    Cell volumes suitability. Donor suitability function corresponds to cell
    volumes in master and donor regions.

    Minimises the overlap by shifting it towards cells with similar size.

SourceFiles
    cellVolumes.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellVolumes_H
#define cellVolumes_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                    Class cellVolumes Declaration
\*---------------------------------------------------------------------------*/

class cellVolumes
:
    public donorSuitability
{

public:

    //- Runtime type information
    TypeName("cellVolumes");


    // Constructors

        //- Construct from components
        cellVolumes
        (
            const overlapFringe& overlapFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellVolumes()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellVolumes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
