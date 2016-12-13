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
    Foam::tetCell

Description
    A tetrahedral cell primitive.

    It is important that the ordering of edges is the same for a tetrahedron
    class, a tetrahedron cell shape model and a tetCell

SourceFiles
    tetCell.C
    tetCellI.H

\*---------------------------------------------------------------------------*/

#ifndef tetCell_H
#define tetCell_H

#include "FixedList.H"
#include "triFace.H"
#include "edge.H"
#include "pointField.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellShape;

/*---------------------------------------------------------------------------*\
                           class tetCell Declaration
\*---------------------------------------------------------------------------*/

class tetCell
:
    public FixedList<label, 4>
{

public:

    // Constructors

        //- Construct null
        inline tetCell();

        //- Construct from four points
        inline tetCell
        (
            const label a,
            const label b,
            const label c,
            const label d
        );

        //- Construct from FixedList
        explicit inline tetCell(const FixedList<label, 4>&);

        //- Construct from Istream
        explicit inline tetCell(Istream&);


    // Member Functions

        // Access

            //- Return i-th face
            inline triFace face(const label faceI) const;

            //- Return first face adjacent to the given edge
            inline label edgeFace(const label edgeI) const;

            //- Return face adjacent to the given face sharing the same edge
            inline label edgeAdjacentFace
            (
                const label edgeI,
                const label faceI
            ) const;

            //- Return i-th edge
            inline edge tetEdge(const label edgeI) const;


        // Operations

            //- Return tet shape cell
            cellShape tetCellShape() const;

            //- Return the tetrahedron
            inline tetPointRef tet(const pointField&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
