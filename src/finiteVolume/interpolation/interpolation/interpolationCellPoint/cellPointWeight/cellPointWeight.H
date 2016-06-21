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
    Foam::cellPointWeight

Description
    Foam::cellPointWeight

SourceFiles
    cellPointWeight.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointWeight_H
#define cellPointWeight_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellPointWeight Declaration
\*---------------------------------------------------------------------------*/

class cellPointWeight
{
protected:

    // Protected data

       //- Cell index
       const label cellIndex_;

       //- Weights applied to tet vertices
       FixedList<scalar, 4> weights_;

       //- Face vertex indices
       FixedList<label, 3> faceVertices_;


    // Protected Member Functions

        void findTetrahedron
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellIndex
        );

        void findTriangle
        (
            const polyMesh& mesh,
            const vector& position,
            const label faceIndex
        );


public:

    //- Debug switch
    static debug::debugSwitch debug;

    //- Tolerance used in calculating barycentric co-ordinates
    //  (applied to normailised values)
    static scalar tol;


    // Constructors

        //- Construct from components
        cellPointWeight
        (
            const polyMesh& mesh,
            const vector& position,
            const label nCell,
            const label facei = -1
        );


    // Member Functions

        //- Cell index
        inline label cell() const
        {
            return cellIndex_;
        }

        //- interpolation weights
        inline const FixedList<scalar, 4>& weights() const
        {
            return weights_;
        }

        //- interpolation addressing for points on face
        inline const FixedList<label, 3>& faceVertices() const
        {
            return faceVertices_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
