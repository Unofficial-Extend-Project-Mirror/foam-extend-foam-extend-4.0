/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "processorSAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorSAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        SAMGInterface,
        processorSAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorSAMGInterface::processorSAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const crMatrix& interfaceProlongation,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const crMatrix& nbrInterfaceProlongation
)
:
    SAMGInterface(lduMesh, interfaceProlongation, nbrInterfaceProlongation),
    fineProcInterface_(refCast<const processorLduInterface>(fineInterface)),
    comm_(fineProcInterface_.comm()),
    tag_(fineProcInterface_.tag())
{
    // READ FIRST
    // This code is written for the FILTERED local matrix
    // because the addressing is very natural. If this is not ok, filter the
    // big matrix here using fineInterface.faceCells(), transpose, and use
    // the same names for crMatrix arrays (then, there is no need to change
    // the triple product code).

    // Algorithm details:
    // Go through the triple product similar to SAMG Policy - collect
    // coarse addressing on master and slave in the same order

    // The resulting contributions from triple product will be saved using
    // the following arrays:
    // - faceCells_: for coarse matrix, it tells us the owner of the
    //               boundary coefficient (on local side)
    // - fineAddressing_: saves the index of the fine boundary coefficient
    //                    in the natural order of appearance on master in
    //                    triple product
    // - restrictAddressing_: saves the coarse face for which the
    //                        contribution (coefficient) needs to sum-up
    // - restrictWeights_: weights for each fine boundary coefficient in
    //                     natural order of appearance on master

    // *Note: (A*B)^T = B^T*A^T, which is exactly what I have. Use the fact
    // on the slave processor!!!

    // Pair collection is always done looking for the master side
    // HJ, 8/Jun/2017

    // Look up neighbours for a master
    HashTable<DynamicList<label, 4>, label, Hash<label> > neighboursTable
    (
        Foam::max(128, fineProcInterface_.interfaceSize())
    );

    // Neighbour face-faces addressing for a face with multiple neighbours
    HashTable<DynamicList<DynamicList<label, 4>, 4>, label, Hash<label> >
    faceFaceTable
    (
        Foam::max(128, fineProcInterface_.interfaceSize())
    );

    // Neighbour face-faces weights for a face with split neighbours
    HashTable<DynamicList<DynamicList<scalar, 4>, 4>, label, Hash<label> >
    faceFaceWeightsTable
    (
        Foam::max(128, fineProcInterface_.interfaceSize())
    );

    // Count the number of coarse faces
    label nCoarseFaces = 0;

    // Count the number of agglomeration pairs
    label nAgglomPairs = 0;

    // Switching prolongation matrices
    const crMatrix* masterP = nullptr;
    const crMatrix* neighbourP = nullptr;

    if (myProcNo() < neighbProcNo())
    {
        // Grab prolongation matrix, master side
        masterP = &interfaceProlongation;
        neighbourP = &nbrInterfaceProlongation;
    }
    else
    {
        // Grab prolongation matrix, slave side
        masterP = &nbrInterfaceProlongation;
        neighbourP = &interfaceProlongation;
    }

    label curMaster, curSlave;

    // Loop through the interface and calculate triple product row-by-row
    for(label faceI = 0; faceI < fineProcInterface_.interfaceSize(); faceI++)
    {
        // Triple product for only one row of prolongation and restriction
        for
        (
            label indexR = masterP->crAddr().rowStart()[faceI];
            indexR < masterP->crAddr().rowStart()[faceI + 1];
            indexR++
        )
        {
            // Grab weight from restriction
            scalar rWeight = masterP->coeffs()[indexR];

            for
            (
                label indexP = neighbourP->crAddr().rowStart()[faceI];
                indexP < neighbourP->crAddr().rowStart()[faceI + 1];
                indexP++
            )
            {
                // Grab weight from prolongation
                scalar pWeight = neighbourP->coeffs()[indexP];

                // Code in the current master and slave - used for
                // identifying the face
                curMaster = masterP->crAddr().column()[indexR];
                curSlave = neighbourP->crAddr().column()[indexP];

                if (neighboursTable.found(curMaster))
                {
                    // This contribution already exists - add the result
                    // to the existing contribution

                    // This master side face already exists

                    // Check all current neighbours to see if the current
                    // slave already exists.  If so, add the coefficient.

                    DynamicList<label, 4>& curNbrs =
                        neighboursTable.find(curMaster)();

                    DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
                        faceFaceTable.find(curMaster)();

                    DynamicList<DynamicList<scalar, 4>, 4>& curFaceWeights =
                        faceFaceWeightsTable.find(curMaster)();

                    // Search for coded neighbour
                    bool nbrFound = false;

                    forAll (curNbrs, curNbrI)
                    {
                        // Check neighbour slave
                        if (curNbrs[curNbrI] == curSlave)
                        {
                            nbrFound = true;
                            curFaceFaces[curNbrI].append(faceI);
                            curFaceWeights[curNbrI].append(pWeight*rWeight);

                            // New agglomeration pair found in already
                            // existing pair
                            nAgglomPairs++;

                            break;
                        }
                    }

                    if (!nbrFound)
                    {
                        curNbrs.append(curSlave);

                        DynamicList<label, 4> newFF;
                        newFF.append(faceI);
                        curFaceFaces.append(newFF);

                        DynamicList<scalar, 4> newFW;
                        newFW.append(pWeight*rWeight);
                        curFaceWeights.append(newFW);

                        // New coarse face created for existing master
                        nCoarseFaces++;
                        nAgglomPairs++;
                    }
                }
                else
                {
                    // This master has got no neighbours yet.
                    // Add a neighbour and a coefficient as a
                    // new list, thus creating a new face

                    DynamicList<label, 4> newNbrs;
                    newNbrs.append(curSlave);
                    neighboursTable.insert
                    (
                        curMaster,
                        newNbrs
                    );

                    DynamicList<DynamicList<label, 4>, 4> newFF;
                    newFF.append(DynamicList<label, 4>());
                    newFF[0].append(faceI);
                    faceFaceTable.insert
                    (
                        curMaster,
                        newFF
                    );

                    DynamicList<DynamicList<scalar, 4>, 4> newFFWeights;
                    newFFWeights.append(DynamicList<scalar, 4>());
                    newFFWeights[0].append(pWeight*rWeight);
                    faceFaceWeightsTable.insert
                    (
                        curMaster,
                        newFFWeights
                    );

                    // New coarse face created for a new master
                    nCoarseFaces++;
                    nAgglomPairs++;
                }
            }
        }
    } // end for all fine faces

    // Since only local faces are analysed, lists can now be resized
    faceCells_.setSize(nCoarseFaces);
    fineAddressing_.setSize(nAgglomPairs);
    restrictAddressing_.setSize(nAgglomPairs);
    restrictWeights_.setSize(nAgglomPairs);

    // Global faces shall be assembled by the increasing label of master
    // cluster ID.
    labelList contents = neighboursTable.toc();

    // Sort makes sure the order is identical on both sides.
    // HJ, 20/Feb/2009 and 6/Jun/2016
    sort(contents);

    // Note: Restriction is done on master side only

    // Reset face counter for re-use
    nCoarseFaces = 0;
    nAgglomPairs = 0;

    // Master side
    if (myProcNo() < neighbProcNo())
    {
        // On master side, the owner addressing is stored in table of contents
        forAll (contents, masterI)
        {
            DynamicList<label, 4>& curNbrs =
                neighboursTable.find(contents[masterI])();

            DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            DynamicList<DynamicList<scalar, 4>, 4>& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            forAll (curNbrs, curNbrI)
            {
                // Get faces and weights
                DynamicList<label, 4>& facesIter = curFaceFaces[curNbrI];

                DynamicList<scalar, 4>& weightsIter = curFaceWeights[curNbrI];

                // Record master cluster index
                faceCells_[nCoarseFaces] = contents[masterI];

                // Collect agglomeration data
                forAll (facesIter, facesIterI)
                {
                    fineAddressing_[nAgglomPairs] = facesIter[facesIterI];
                    restrictAddressing_[nAgglomPairs] = nCoarseFaces;
                    restrictWeights_[nAgglomPairs] = weightsIter[facesIterI];

                    nAgglomPairs++;
                }

                nCoarseFaces++;
            }
        }
    }
    // Slave side
    else
    {
        // On slave side, the owner addressing is stored in linked lists
        forAll (contents, masterI)
        {
            DynamicList<label, 4>& curNbrs =
                neighboursTable.find(contents[masterI])();

            DynamicList<DynamicList<label, 4>, 4>& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            DynamicList<DynamicList<scalar, 4>, 4>& curFaceWeights =
                faceFaceWeightsTable.find(contents[masterI])();

            forAll (curNbrs, curNbrI)
            {
                // Get faces and weights
                DynamicList<label, 4>& facesIter = curFaceFaces[curNbrI];

                DynamicList<scalar, 4>& weightsIter = curFaceWeights[curNbrI];

                faceCells_[nCoarseFaces] = curNbrs[curNbrI];

                forAll (facesIter, facesIterI)
                {
                    fineAddressing_[nAgglomPairs] = facesIter[facesIterI];
                    restrictAddressing_[nAgglomPairs] = nCoarseFaces;
                    restrictWeights_[nAgglomPairs] = weightsIter[facesIterI];

                    nAgglomPairs++;
                }

                nCoarseFaces++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::processorSAMGInterface::~processorSAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorSAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    send(commsType, interfaceData);
}


Foam::tmp<Foam::labelField> Foam::processorSAMGInterface::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorSAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    send(commsType, interfaceInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorSAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorSAMGInterface::initProlongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    // Send prolongation matrix, using IOstream operators
    OPstream toNbr(Pstream::blocking, neighbProcNo());
    toNbr<< filteredP;
}


Foam::autoPtr<Foam::crMatrix>
Foam::processorSAMGInterface::prolongationTransfer
(
    const Pstream::commsTypes commsType,
    const crMatrix& filteredP
) const
{
    IPstream fromNbr(Pstream::blocking, neighbProcNo());

    autoPtr<crMatrix> tnbrP(new crMatrix(fromNbr));

    return tnbrP;
}


// ************************************************************************* //
