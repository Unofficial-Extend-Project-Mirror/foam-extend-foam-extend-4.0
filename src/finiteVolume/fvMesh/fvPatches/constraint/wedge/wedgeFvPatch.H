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
    Foam::wedgeFvPatch

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFvPatch_H
#define wedgeFvPatch_H

#include "fvPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wedgeFvPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeFvPatch
:
    public fvPatch
{
    // Private data

        const wedgePolyPatch& wedgePolyPatch_;

protected:

    // Protected Member Functions

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const;


public:

    //- Runtime type information
    TypeName(wedgePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        wedgeFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm);


    // Destructor

        virtual ~wedgeFvPatch()
        {}


    // Member functions

        // Access

            //- Return face transformation tensor
            const tensor& faceT() const
            {
                return wedgePolyPatch_.faceT();
            }

            //- Return neighbour-cell transformation tensor
            const tensor& cellT() const
            {
                return wedgePolyPatch_.cellT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
