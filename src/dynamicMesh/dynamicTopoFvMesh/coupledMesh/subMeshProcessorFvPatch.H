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
    subMeshProcessorFvPatch

Description
    Customized processor type for subMeshes.

    The intent of this class is to provide processor patch information
    on subMeshes, but avoiding all forms of communication.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef subMeshProcessorFvPatch_H
#define subMeshProcessorFvPatch_H

#include "processorFvPatch.H"
#include "subMeshProcessorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class subMeshProcessorFvPatch Declaration
\*---------------------------------------------------------------------------*/

class subMeshProcessorFvPatch
:
    public processorFvPatch
{
    // Private Data

        const subMeshProcessorPolyPatch& procPolyPatch_;

protected:

    // Protected Member Functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;

public:

    //- Runtime type information
    TypeName(subMeshProcessorPolyPatch::typeName_());

    // Constructors

        //- Construct from components
        subMeshProcessorFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            processorFvPatch(patch, bm),
            procPolyPatch_(refCast<const subMeshProcessorPolyPatch>(patch))
        {}

    // Destructor

        virtual ~subMeshProcessorFvPatch()
        {}

    // Member functions

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
