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
    tetFemSolution

Description
    Selector class for finite element solution solution.
    tetFemMesh is derived from tetFemSolution so that all fields
    have access to the tetFemSolution from the mesh reference they
    hold.

\*---------------------------------------------------------------------------*/

#ifndef tetFemSolution_H
#define tetFemSolution_H

#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tetFemSolution Declaration
\*---------------------------------------------------------------------------*/

class tetFemSolution
:
    public solution
{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        tetFemSolution(const tetFemSolution&);
        void operator=(const tetFemSolution&);


public:

    // Constructors

        //- Construct from objectRegistry
        tetFemSolution(const objectRegistry& obr)
        :
            solution(obr, "tetFemSolution")
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
