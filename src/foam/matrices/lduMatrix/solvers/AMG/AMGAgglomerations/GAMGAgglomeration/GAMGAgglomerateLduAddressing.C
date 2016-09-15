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

\*---------------------------------------------------------------------------*/

#include "GAMGAgglomeration.H"
#include "AMGInterface.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GAMGAgglomeration::agglomerateLduAddressing
(
    const label fineLevelIndex
)
{
    const lduMesh& fineMesh = meshLevel(fineLevelIndex);
    const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

    const unallocLabelList& upperAddr = fineMeshAddr.upperAddr();
    const unallocLabelList& lowerAddr = fineMeshAddr.lowerAddr();

    label nFineFaces = upperAddr.size();

    // Get restriction map for current level
    const labelField& restrictMap = restrictAddressing(fineLevelIndex);

    if (min(restrictMap) == -1)
    {
        FatalErrorIn("GAMGAgglomeration::agglomerateLduAddressing")
            << "min(restrictMap) == -1" << exit(FatalError);
    }

    if (restrictMap.size() != fineMeshAddr.size())
    {
        FatalErrorIn
        (
            "GAMGAgglomeration::agglomerateLduAddressing"
            "(const label fineLevelIndex)"
        )   << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << restrictMap.size()
            << " nEqns: " << fineMeshAddr.size()
            << abort(FatalError);
    }


    // Get the number of coarse cells
    const label nCoarseCells = nCells_[fineLevelIndex];

    // Storage for coarse cell neighbours and coefficients

    // Guess initial maximum number of neighbours in coarse cell
    label maxNnbrs = 10;

    // Number of faces for each coarse-cell
    labelList cCellnFaces(nCoarseCells, 0);

    // Setup initial packed storage for coarse-cell faces
    labelList cCellFaces(maxNnbrs*nCoarseCells);

    // Create face-restriction addressing
    faceRestrictAddressing_.set(fineLevelIndex, new labelList(nFineFaces));
    labelList& faceRestrictAddr = faceRestrictAddressing_[fineLevelIndex];

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineFaces);

    // Counter for coarse faces
    label nCoarseFaces = 0;

    // Loop through all fine faces
    forAll (upperAddr, fineFacei)
    {
        label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
        label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine face inside of a coarse cell keep the address
            // of the cell corresponding to the face in the faceRestrictAddr
            // as a negative index
            faceRestrictAddr[fineFacei] = -(rmUpperAddr + 1);
        }
        else
        {
            // this face is a part of a coarse face

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // check the neighbour to see if this face has already been found
            label* ccFaces = &cCellFaces[maxNnbrs*cOwn];

            bool nbrFound = false;
            label& ccnFaces = cCellnFaces[cOwn];

            for (int i=0; i<ccnFaces; i++)
            {
                if (initCoarseNeighb[ccFaces[i]] == cNei)
                {
                    nbrFound = true;
                    faceRestrictAddr[fineFacei] = ccFaces[i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnFaces >= maxNnbrs)
                {
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    cCellFaces.setSize(maxNnbrs*nCoarseCells);

                    forAllReverse(cCellnFaces, i)
                    {
                        label* oldCcNbrs = &cCellFaces[oldMaxNnbrs*i];
                        label* newCcNbrs = &cCellFaces[maxNnbrs*i];

                        for (int j=0; j<cCellnFaces[i]; j++)
                        {
                            newCcNbrs[j] = oldCcNbrs[j];
                        }
                    }

                    ccFaces = &cCellFaces[maxNnbrs*cOwn];
                }

                ccFaces[ccnFaces] = nCoarseFaces;
                initCoarseNeighb[nCoarseFaces] = cNei;
                faceRestrictAddr[fineFacei] = nCoarseFaces;
                ccnFaces++;

                // new coarse face created
                nCoarseFaces++;
            }
        }
    } // end for all fine faces


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    labelList coarseOwner(nCoarseFaces);
    labelList coarseNeighbour(nCoarseFaces);
    labelList coarseFaceMap(nCoarseFaces);

    label coarseFacei = 0;

    forAll (cCellnFaces, cci)
    {
        label* cFaces = &cCellFaces[maxNnbrs*cci];
        label ccnFaces = cCellnFaces[cci];

        for (int i=0; i<ccnFaces; i++)
        {
            coarseOwner[coarseFacei] = cci;
            coarseNeighbour[coarseFacei] = initCoarseNeighb[cFaces[i]];
            coarseFaceMap[cFaces[i]] = coarseFacei;
            coarseFacei++;
        }
    }

    forAll(faceRestrictAddr, fineFacei)
    {
        if (faceRestrictAddr[fineFacei] >= 0)
        {
            faceRestrictAddr[fineFacei] =
                coarseFaceMap[faceRestrictAddr[fineFacei]];
        }
    }

    // Clear the temporary storage for the coarse cell data
    cCellnFaces.setSize(0);
    cCellFaces.setSize(0);
    initCoarseNeighb.setSize(0);
    coarseFaceMap.setSize(0);


    // Create coarse-level interfaces

    // Get reference to fine-level interfaces
    const lduInterfacePtrsList& fineInterfaces =
        interfaceLevels_[fineLevelIndex];

    // Initialise transfer of restrict addressing on the interface
    forAll (fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            if (fineInterfaces[inti].coupled())
            {
                fineInterfaces[inti].initInternalFieldTransfer
                (
                    Pstream::blocking,
                    restrictMap
                );
            }
        }
    }

    // Store coefficients to avoid tangled communications
    // HJ, 1/Apr/2009
    FieldField<Field, label> fineInterfaceAddr(fineInterfaces.size());

    forAll (fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            if (fineInterfaces[inti].coupled())
            {
                fineInterfaceAddr.set
                (
                    inti,
                    new labelField
                    (
                        fineInterfaces[inti].internalFieldTransfer
                        (
                            Pstream::blocking,
                            restrictMap
                        )
                    )
                );
            }
        }
    }

    // Add the coarse level

    // Set the coarse ldu addressing onto the list
    meshLevels_.set
    (
        fineLevelIndex,
        new lduPrimitiveMesh
        (
            nCoarseCells,
            coarseOwner,
            coarseNeighbour,
            fineMesh.comm(),
            true
        )
    );

    // Create coarse-level interfaces
    interfaceLevels_.set
    (
        fineLevelIndex + 1,
        new lduInterfacePtrsList(fineInterfaces.size())
    );

    lduInterfacePtrsList& coarseInterfaces =
        interfaceLevels_[fineLevelIndex + 1];

    labelListList coarseInterfaceAddr(fineInterfaces.size());

    forAll (fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            if (fineInterfaces[inti].coupled())
            {
                coarseInterfaces.set
                (
                    inti,
                    AMGInterface::New
                    (
                        meshLevels_[fineLevelIndex],
                        coarseInterfaces,
                        fineInterfaces[inti],
                        fineInterfaces[inti].interfaceInternalField
                        (
                            restrictMap
                        ),
                        fineInterfaceAddr[inti]
                    ).ptr()
                );

                coarseInterfaceAddr[inti] = coarseInterfaces[inti].faceCells();
            }
        }
    }

    meshLevels_[fineLevelIndex].addInterfaces
    (
        coarseInterfaces,
        coarseInterfaceAddr,
        fineMeshAddr.patchSchedule()
    );
}


// ************************************************************************* //
