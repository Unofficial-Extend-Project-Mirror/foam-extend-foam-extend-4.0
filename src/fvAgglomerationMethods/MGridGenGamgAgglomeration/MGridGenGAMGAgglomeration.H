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
    Foam::MGridGenGAMGAgglomeration

Description
    Agglomerate using the MGridGen algorithm.

SourceFiles
    MGridGenGAMGAgglomeration.C
    MGridGenGAMGAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef MGridGenGAMGAgglomeration_H
#define MGridGenGAMGAgglomeration_H

#include "fvMesh.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class MGridGenGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class MGridGenGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        const fvMesh& fvMesh_;


    // Private Member Functions

        //- Construct the CSR format addressing
        void makeCompactCellFaceAddressingAndFaceWeights
        (
            const lduAddressing& fineAddressing,
            labelList& cellCells,
            labelList& cellCellOffsets,
            const vectorField& Si,
            List<scalar>& faceWeights
        );

        //- Calculate and return agglomeration
        tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const label minSize,
            const label maxSize,
            const lduAddressing& fineAddressing,
            const scalarField& V,
            const vectorField& Sf,
            const scalarField& Sb
        );


        //- Disallow default bitwise copy construct
        MGridGenGAMGAgglomeration(const MGridGenGAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const MGridGenGAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("MGridGen");


    // Constructors

        //- Construct given mesh and controls
        MGridGenGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
