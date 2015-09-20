/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
{
    // Private data

        //- Solver settings
        solution solutionDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        tetFemSolution(const tetFemSolution&);
        void operator=(const tetFemSolution&);


public:

    // Constructors

        //- Construct from objectRegistry
        tetFemSolution(const objectRegistry& obr)
        :
            solutionDict_(obr, "tetFemSolution")
        {}


    // Member functions

        //- Return access to solver settings
        const solution& solutionDict() const
        {
            return solutionDict_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
