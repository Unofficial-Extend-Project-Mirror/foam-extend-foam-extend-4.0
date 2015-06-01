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
    Foam::symmetryFaPatch

Description
    Symmetry-plane patch.

SourceFiles
    symmetryFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryFaPatch_H
#define symmetryFaPatch_H

#include "faPatch.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class symmetryFaPatch Declaration
\*---------------------------------------------------------------------------*/

class symmetryFaPatch
:
    public faPatch
{

protected:

    // Protected Member Functions

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const;


public:

    //- Runtime type information
    TypeName(symmetryPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        symmetryFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex
        );

        //- Construct from dictionary
        symmetryFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );

        //- Construct and return a clone, resetting the edge list
        //  and boundary mesh
        virtual autoPtr<faPatch> clone
        (
            const faBoundaryMesh& bm,
            const labelList& edgeLabels,
            const label index,
            const label ngbPolyPatchIndex
        ) const
        {
            return autoPtr<faPatch>
            (
                new symmetryFaPatch
                (
                    name(), edgeLabels, index, bm, ngbPolyPatchIndex
                )
            );
        }

    // Destructor

        virtual ~symmetryFaPatch()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
