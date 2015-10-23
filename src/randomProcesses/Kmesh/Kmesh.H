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
    Foam::Kmesh

Description
    Calculate the wavenumber vector field corresponding to the
    space vector field of a finite volume mesh;

SourceFiles
    Kmesh.C

\*---------------------------------------------------------------------------*/

#ifndef Kmesh_H
#define Kmesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Kmesh Declaration
\*---------------------------------------------------------------------------*/

class Kmesh
:
    public vectorField
{

    // Private data

        //- Dimensions of box
        vector L;

        //- Multi-dimensional addressing array
        labelList NN;

        //- Maximum wavenumber
        scalar Kmax;

public:

    // Constructors

        //- Construct from fvMesh
        Kmesh(const fvMesh&);


    // Member Functions

        // Access

        const vector& sizeOfBox() const
        {
            return L;
        }

        const labelList& nn() const
        {
            return NN;
        }

        scalar max() const
        {
            return Kmax;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
