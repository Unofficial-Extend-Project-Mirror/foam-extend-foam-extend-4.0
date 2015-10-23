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
    Foam::block

Description
    Creates a single block of cells from point coordinates,
    numbers of cells in each direction and expansion ratio


SourceFiles
    block.C
    blockIO.C


\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "pointField.H"
#include "labelList.H"

#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
{
    // Private data

        //- block definition
        blockDescriptor BlockDef;

        //- list of vertices
        pointField Vertices;

        //- list of cells
        labelListList Cells;

        //- boundary patches
        labelListListList BoundaryPatches;


    // private member functions

        label vtxLabel(label i, label j, label k);

        void blockPoints();

        void blockCells();

        void blockBoundary();

public:

    // Constructors

        //- from the block definition
        block(const blockDescriptor&);

        //- as copy
        block(const block&);

        //- clone function
        autoPtr<block> clone() const
        {
            notImplemented("block::clone()");
            return autoPtr<block>(NULL);
        }

    // Destructor

        ~block();


    // Member Functions

        // Access

        const blockDescriptor& blockDef() const;
        const pointField& points() const;
        const labelListList& cells() const;
        const labelListListList& boundaryPatches() const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const block&);
};


inline Istream& operator>>(Istream& is, block*)
{
    notImplemented("Istream& operator>>(Istream& is, block*)");
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
