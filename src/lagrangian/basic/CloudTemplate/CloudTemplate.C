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

#include "CloudTemplate.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "foamTime.H"
#include "OFstream.H"

#include "profiling.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    particleCount_(0)
{
    IDLList<ParticleType>::operator=(particles);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    particleCount_(0)
{
    IDLList<ParticleType>::operator=(particles);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::label Foam::Cloud<ParticleType>::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningIn("Cloud<ParticleType>::getNewParticleID() const")
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


namespace Foam
{

class combineNsTransPs
{

public:

    void operator()(labelListList& x, const labelListList& y) const
    {
        forAll(y, i)
        {
            if (y[i].size())
            {
                x[i] = y[i];
            }
        }
    }
};

} // End namespace Foam


template<class ParticleType>
template<class TrackingData>
void Foam::Cloud<ParticleType>::move(TrackingData& td)
{
    addProfile2(moveProfile,"Cloud<ParticleType>::move_"+this->name());

    const globalMeshData& pData = polyMesh_.globalData();
    const labelList& processorPatches = pData.processorPatches();
    const labelList& processorPatchIndices = pData.processorPatchIndices();
    const labelList& processorPatchNeighbours =
        pData.processorPatchNeighbours();

    // Initialise the setpFraction moved for the particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().stepFraction() = 0;
    }

    // Assume there will be particles to transfer
    bool transfered = true;

    // While there are particles to transfer
    while (transfered)
    {
        // List of lists of particles to be transfered for all the processor
        // patches
        List<IDLList<ParticleType> > transferList(processorPatches.size());

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            bool keepParticle = p.move(td);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if (Pstream::parRun() && p.facei_ >= pMesh().nInternalFaces())
                {
                    label patchi = pMesh().boundaryMesh().whichPatch(p.facei_);
                    label n = processorPatchIndices[patchi];

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (n != -1)
                    {
                        p.prepareForParallelTransfer(patchi, td);
                        transferList[n].append(this->remove(&p));
                    }
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (Pstream::parRun())
        {
            // List of the numbers of particles to be transfered across the
            // processor patches
            labelList nsTransPs(transferList.size());

            forAll(transferList, i)
            {
                nsTransPs[i] = transferList[i].size();
            }

            // List of the numbers of particles to be transfered across the
            // processor patches for all the processors
            labelListList allNTrans(Pstream::nProcs());
            allNTrans[Pstream::myProcNo()] = nsTransPs;
            combineReduce(allNTrans, combineNsTransPs());

            transfered = false;

            forAll(allNTrans, i)
            {
                forAll(allNTrans[i], j)
                {
                    if (allNTrans[i][j])
                    {
                        transfered = true;
                        break;
                    }
                }
            }

            if (!transfered)
            {
                break;
            }

            forAll(transferList, i)
            {
                if (transferList[i].size())
                {
                    OPstream particleStream
                    (
                        Pstream::blocking,
                        refCast<const processorPolyPatch>
                        (
                            pMesh().boundaryMesh()[processorPatches[i]]
                        ).neighbProcNo()
                    );

                    particleStream << transferList[i];
                }
            }

            forAll(processorPatches, i)
            {
                label patchi = processorPatches[i];

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (pMesh().boundaryMesh()[patchi]);

                label neighbProci =
                    procPatch.neighbProcNo() - Pstream::masterNo();

                label neighbProcPatchi = processorPatchNeighbours[patchi];

                label nRecPs = allNTrans[neighbProci][neighbProcPatchi];

                if (nRecPs)
                {
                    IPstream particleStream
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    IDLList<ParticleType> newParticles
                    (
                        particleStream,
                        typename ParticleType::iNew(*this)
                    );

                    forAllIter
                    (
                        typename Cloud<ParticleType>,
                        newParticles,
                        newpIter
                    )
                    {
                        ParticleType& newp = newpIter();
                        newp.correctAfterParallelTransfer(patchi, td);
                        addParticle(newParticles.remove(&newp));
                    }
                }
            }
        }
        else
        {
            transfered = false;
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (cloud::debug)
    {
        Info<< "Cloud<ParticleType>::autoMap(const morphFieldMapper& map) "
               "for lagrangian cloud " << cloud::name() << endl;
    }

    const labelList& reverseCellMap = mapper.reverseCellMap();
    const labelList& reverseFaceMap = mapper.reverseFaceMap();

    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        if (reverseCellMap[pIter().celli_] >= 0)
        {
            pIter().celli_ = reverseCellMap[pIter().celli_];

            if (pIter().facei_ >= 0 && reverseFaceMap[pIter().facei_] >= 0)
            {
                pIter().facei_ = reverseFaceMap[pIter().facei_];
            }
            else
            {
                pIter().facei_ = -1;
            }
        }
        else
        {
            label trackStartCell = mapper.mergedCell(pIter().celli_);

//          Tommaso, bug fix 27/02/2008
//          when the topology of the mesh is changed
//          it is necessary to update the cell in which the parcel is,
//          in particular when cells are removed

            // HJ, merge 1.5.x  20/Oct/2008
            if (trackStartCell < 0)
            {
                pIter().celli_ = polyMesh_.findCell(pIter().position());
            }

            vector p = pIter().position();
            const_cast<vector&>(pIter().position()) =
                polyMesh_.cellCentres()[trackStartCell];
            pIter().stepFraction() = 0;
            pIter().track(p);
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::split
(
    const labelList& cellToProc,
    const labelListList& procCellAddressing,
    const labelListList& procFaceAddressing,
    List<IDLList<ParticleType> >& procParticles
)
{
    // Create addressing old cell ID to new cell ID from procCellAddressing
    labelList reverseCellMap(cellToProc.size());
    forAll(procCellAddressing, procI)
    {
        const labelList& procCell = procCellAddressing[procI];

        forAll(procCell, cellI)
        {
            reverseCellMap[procCell[cellI]] = cellI;
        }
    }

    labelList reverseFaceMap(polyMesh_.nFaces());
    forAll(procFaceAddressing, procI)
    {
        const labelList& procFace = procFaceAddressing[procI];

        forAll(procFace, faceI)
        {
            reverseFaceMap[mag(procFace[faceI])-1] = faceI;
        }
    }

    if (cloud::debug)
    {
        Pout << "printing cloud before splitting" << endl;

        forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
        {
            const ParticleType& p = pIter();

            Pout << p.celli_ << " "
                << p.facei_ << " "
                << p.position() << " "
                << p.inCell() << endl;
        }
    }

    // Loop over all particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();

        const label newProc = cellToProc[p.celli_];

        // Map cellID and faceID of the particle to new mesh
        p.celli_ = reverseCellMap[p.celli_];

        if (p.facei_ >= 0)
        {
            // HR 29.06.17: face mapping not tested. Did not occur in test cases.
            p.facei_ = reverseFaceMap[p.facei_];
        }

        if (newProc != Pstream::myProcNo())
        {
            // Add for transfer to new processor and delete from local cloud
            procParticles[newProc].append(this->remove(&p));
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::rebuild
(
    List<IDLList<ParticleType> >& receivedClouds,
    const PtrList<labelIOList>& cellProcAddressing,
    const PtrList<labelIOList>& faceProcAddressing
)
{
    if (cellProcAddressing.set(Pstream::myProcNo()))
    {
        const labelList& procCell = cellProcAddressing[Pstream::myProcNo()];
        const labelList& procFace = faceProcAddressing[Pstream::myProcNo()];

        // Particles in local cloud. Map cellID and faceID in-place
        forAllIter
        (
            typename Cloud<ParticleType>,
            *this,
            pIter
        )
        {
            ParticleType& p = pIter();

            // Map cellID and faceID of the particle to new mesh
            p.celli_ = procCell[p.celli_];

            if (p.facei_ >= 0)
            {
                // HR 29.06.17: face mapping not tested. Did not occur in test cases.
                p.facei_ = procFace[p.facei_];
            }
        }
    }

    forAll(cellProcAddressing, procI)
    {
        if (!cellProcAddressing.set(procI))
        {
            continue;
        }

        const labelList& procCell = cellProcAddressing[procI];
        const labelList& procFace = faceProcAddressing[procI];

        IDLList<ParticleType>& receivedCloud = receivedClouds[procI];

        // Particles received from other processors. To be added to local cloud
        forAllIter
        (
            typename Cloud<ParticleType>,
            receivedCloud,
            newpIter
        )
        {
            ParticleType& newp = newpIter();

            // Map cellID and faceID of the particle to new mesh
            newp.celli_ = procCell[newp.celli_];

            if (newp.facei_ >= 0)
            {
                // HR 29.06.17: face mapping not tested. Did not occur in test cases.
                newp.facei_ = procFace[newp.facei_];
            }

            // Add to cloud and remove from transfer list
            this->addParticle(receivedCloud.remove(&newp));
        }
    }


    if (cloud::debug)
    {
        Pout << "printing cloud after rebuild" << endl;

        forAllIter
        (
            typename Cloud<ParticleType>,
            *this,
            pIter
        )
        {
            const ParticleType& p = pIter();

            Pout << p.celli_ << " "
                << p.facei_ << " "
                << p.position() << " "
                << p.inCell() << endl;
        }
    }

    this->updateMesh();
}


template<class ParticleType>
Foam::labelList Foam::Cloud<ParticleType>::nParticlesPerCell() const
{
    Pout << "In nParticlesPerCell" << endl;
    labelList nppc (pMesh().nCells(), 0);

    forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        const label celli = p.cell();

        // Check
        if (celli < 0 || celli >= pMesh().nCells())
        {
            FatalErrorIn
            (
                "Foam::Cloud<ParticleType>::nParticlesPerCell()"
            )
                << "Illegal cell number " << celli
                << " at position " << p.position() << nl
                << "Cell number should be between 0 and "
                << pMesh().nCells()-1 << nl
                << "On this mesh the particle should be in cell "
                << pMesh().findCell(p.position())
                << exit(FatalError);
        }

        nppc[celli]++;
    }

    return nppc;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudTemplateIO.C"

// ************************************************************************* //
