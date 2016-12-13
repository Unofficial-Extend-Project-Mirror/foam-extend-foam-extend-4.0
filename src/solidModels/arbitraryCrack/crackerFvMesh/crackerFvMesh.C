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

#include "crackerFvMesh.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(crackerFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, crackerFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crackerFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

//     if (topoChanger_.size() || faceZones().size() > 0)
//     {
//         Info<< "void crackerFvMesh::addZonesAndModifiers() : "
//             << "Zones and modifiers already present.  Skipping."
//             << endl;

//         return;
//     }

    Info<< "Time = " << time().timeName() << endl
        << "Adding topo modifier to the mesh" << endl;

    const word crackPatchName(dict_.lookup("crackPatch"));
    label crackPatchIndex = boundaryMesh().findPatchID(crackPatchName);

    if (crackPatchIndex < 0)
    {
        FatalErrorIn("void crackerFvMesh::addZonesAndModifiers() const")
            << "Crack patch not found in boundary"
            << abort(FatalError);
    }

    const word openPatchName(dict_.lookup("openPatch"));
    label openPatchIndex = boundaryMesh().findPatchID(openPatchName);

    if (openPatchIndex < 0)
    {
        FatalErrorIn("void crackerFvMesh::addZonesAndModifiers() const")
            << "Open patch not found in boundary"
            << abort(FatalError);
    }

    // Add zones
    if (faceZones().size() == 0)
    {
        List<pointZone*> pz(0);
        List<faceZone*> fz(1);

        fz[0] = new faceZone
        (
            crackPatchName + "Zone",
            labelList(0),
            boolList(0),
            0,
            faceZones()
        );

        List<cellZone*> cz(0);


        Info << "Adding point, face and cell zones" << endl;
        addZones(pz, fz, cz);
    }
    else if (faceZones().findZoneID(crackPatchName + "Zone") == -1)
    {
//         List<pointZone*> pz(0);
//         List<faceZone*> fz(faceZones().size() + 1);

//         fz[0] =
//             new faceZone
//             (
//                 crackPatchName + "Zone",
//                 labelList(0),
//                 boolList(0),
//                 0,
//                 faceZones()
//             );

//         for (label i=0; i<faceZones().size(); i++)
//         {
//             fz[i+1] = new faceZone
//             (
//                 faceZones()[i].name(),
//                 faceZones()[i],
//                 faceZones()[i].flipMap(),
//                 i+1,
//                 faceZones()
//             );
//         }

//         List<cellZone*> cz(0);

//         Info << "Adding crack face zone. "
//             << "Preserving existing face zones." << endl;

//         if (pointZones().size() > 0)
//         {
//             pointZones().clear();
//         }
//         if (faceZones().size() > 0)
//         {
//             faceZones().clear();
//         }
//         if (cellZones().size() > 0)
//         {
//             cellZones().clear();
//         }
//         addZones(pz, fz, cz);

        FatalErrorIn("void crackerFvMesh::addZonesAndModifiers() const")
            << "Crack face zone must be added by hand" << nl
            << " crackZone " << nl
            << "{" << nl
            << "\ttype faceZone;" << nl
            << "\tfaceLabels      0();" << nl
            << "\t  flipMap         0();" << nl
            << "}" << nl
            << abort(FatalError);
    }
    else
    {
        Info << "Face zones already present" << endl;
    }

    // Add a topology modifier
    if (topoChanger_.size() == 0)
    {
        Info << "Adding topology modifiers" << endl;
        topoChanger_.setSize(1);
        topoChanger_.set
        (
            0,
            new faceCracker
            (
                "cracker",
                0,
                topoChanger_,
                crackPatchName + "Zone",
                crackPatchName,
                openPatchName
            )
        );

        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        Info<< "void crackerFvMesh::addZonesAndModifiers() : "
            << "Modifiers already present."
            << endl;
    }

    // Write mesh
    write();
}

void Foam::crackerFvMesh::makeRegions() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (regionsPtr_)
    {
        FatalErrorIn("crackerFvMesh::makeRegions() const")
            << "regions already exist"
            << abort(FatalError);
    }

    regionsPtr_ = new regionSplit(*this);
}

void Foam::crackerFvMesh::makeNCellsInRegion() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (nCellsInRegionPtr_)
    {
        FatalErrorIn("crackerFvMesh::makeNCellsInRegion() const")
            << "number of cells in regions already exist"
            << abort(FatalError);
    }

    nCellsInRegionPtr_ = new labelList(regions().nRegions(), 0);

    labelList& nCellsInRegion = *nCellsInRegionPtr_;

    const labelList& regs = regions();

    forAll(regs, cellI)
    {
        nCellsInRegion[regs[cellI]]++;
    }
}

void Foam::crackerFvMesh::makeGlobalCrackFaceCentresAndSizes() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalCrackFaceCentresPtr_ || globalCrackFaceSizesPtr_)
    {
        FatalErrorIn
        (
            "crackerFvMesh::makeGlobalCrackFaceCentresAndSizes() const"
        )
            << "global crack face centres and sizes already exist"
            << abort(FatalError);
    }


    // Number of faces in global crack
    labelList sizes(Pstream::nProcs(), 0);
    sizes[Pstream::myProcNo()] = boundaryMesh()[crackPatchID_.index()].size();

    if (Pstream::parRun())
    {
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << sizes[Pstream::myProcNo()];
                }
            }
        }

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> sizes[procI];
                }
            }
        }
    }

    label globalCrackSize = sum(sizes);

    globalCrackFaceCentresPtr_ =
        new vectorField(globalCrackSize, vector::zero);
    vectorField& globalCrackFaceCentres = *globalCrackFaceCentresPtr_;

    globalCrackFaceSizesPtr_ =
        new scalarField(globalCrackSize, 0);
    scalarField& globalCrackFaceSizes = *globalCrackFaceSizesPtr_;

    localCrackStart_ = 0;
    for (label procI = 0; procI < Pstream::myProcNo(); procI++)
    {
        localCrackStart_ += sizes[procI];
    }

//     const vectorField& crackCf =
//         boundaryMesh()[crackPatchID_.index()].faceCentres();
    const vectorField::subField crackCf =
        boundaryMesh()[crackPatchID_.index()].faceCentres();

    // Calc face sizes
//     const vectorField& crackSf =
//         boundaryMesh()[crackPatchID_.index()].faceAreas();
    const vectorField::subField crackSf =
        boundaryMesh()[crackPatchID_.index()].faceAreas();

    scalarField delta(crackSf.size(), 0);

    if (nGeometricD() == 3)
    {
        delta = Foam::sqrt(mag(crackSf));
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = bounds().span()[dir];
                break;
            }
        }

        delta = mag(crackSf)/thickness;
    }

    label j=0;
    for
    (
        label i=localCrackStart_;
        i<(localCrackStart_ + sizes[Pstream::myProcNo()]);
        i++
    )
    {
        globalCrackFaceCentres[i] = crackCf[j];
        globalCrackFaceSizes[i] = delta[j];
        j++;
    }

    // Parallel data exchange: collect crack face centres and sizes
    // on all processors
    reduce(globalCrackFaceCentres, sumOp<vectorField>());
    reduce(globalCrackFaceSizes, sumOp<scalarField>());
}


void Foam::crackerFvMesh::makeGlobalCrackFaceAddressing() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (globalCrackFaceAddressingPtr_)
    {
        FatalErrorIn("crackerFvMesh::makeGlobalCrackFaceAddressing() const")
            << "global crack face addressing already exists"
            << abort(FatalError);
    }

    const vectorField& gcfc = globalCrackFaceCentres();
    const scalarField& gcfs = globalCrackFaceSizes();

    globalCrackFaceAddressingPtr_ = new labelList(gcfc.size(), -1);
    labelList& gcfa = *globalCrackFaceAddressingPtr_;

    forAll(gcfa, faceI)
    {
        if (gcfa[faceI] < 0)
        {
            forAll(gcfc, fI)
            {
                if ((fI != faceI) && (gcfa[fI] < 0))
                {
                    if (mag(gcfc[faceI] - gcfc[fI]) < 1e-3*gcfs[faceI])
                    {
                        gcfa[faceI] = fI;
                        gcfa[fI] = faceI;
                        break;
                    }
                }
            }
        }
    }

    // Check addressing
    forAll(gcfa, faceI)
    {
        if (gcfa[faceI] < 0)
        {
            FatalErrorIn
            (
                "crackerFvMesh::makeGlobalCrackFaceAddressing() const"
            )
            << "problem with defining global crack face addressing"
            << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::crackerFvMesh::crackerFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    topoChangeMap_(),
    crackPatchID_(word(dict_.lookup("crackPatch")), boundaryMesh()),
    regionsPtr_(NULL),
    nCellsInRegionPtr_(NULL),
    globalCrackFaceCentresPtr_(NULL),
    globalCrackFaceSizesPtr_(NULL),
    localCrackStart_(-1),
    globalCrackFaceAddressingPtr_(NULL)
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::crackerFvMesh::~crackerFvMesh()
{
    deleteDemandDrivenData(regionsPtr_);
    deleteDemandDrivenData(nCellsInRegionPtr_);
    deleteDemandDrivenData(globalCrackFaceCentresPtr_);
    deleteDemandDrivenData(globalCrackFaceSizesPtr_);
    deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::crackerFvMesh::setBreak
(
    const labelList& facesToBreak,
    const boolList& faceFlip,
    const labelList& coupledFacesToBreak
)
{
    faceCracker& fc = refCast<faceCracker>(topoChanger_[0]);

    fc.setBreak(facesToBreak, faceFlip, coupledFacesToBreak);
}


bool Foam::crackerFvMesh::update()
{
//     autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
    topoChangeMap_ = topoChanger_.changeMesh();

    if (topoChangeMap_.valid())
    {
        deleteDemandDrivenData(regionsPtr_);
        deleteDemandDrivenData(nCellsInRegionPtr_);
        deleteDemandDrivenData(globalCrackFaceCentresPtr_);
        deleteDemandDrivenData(globalCrackFaceSizesPtr_);
        localCrackStart_ = -1;
        deleteDemandDrivenData(globalCrackFaceAddressingPtr_);
    }

//     const labelList& faceMap = topoChangeMap().faceMap();

//     forAll(this->boundaryMesh(), patchI)
//     {
//         label start  = this->boundaryMesh()[patchI].start();
//         label size  = this->boundaryMesh()[patchI].size();

//         Info << this->boundaryMesh()[patchI].name() << endl;
//         for (label i=0;i<size;i++)
//         {
//             Info << faceMap[start+i] << endl;
//         }
//     }

    return topoChangeMap_.valid();
}


const Foam::regionSplit& Foam::crackerFvMesh::regions() const
{
    if (!regionsPtr_)
    {
        makeRegions();
    }

    return *regionsPtr_;
}


Foam::label Foam::crackerFvMesh::nCellsInRegion(label regI) const
{
    if (!nCellsInRegionPtr_)
    {
        makeNCellsInRegion();
    }

    if ((regI < 0) || (regI >= regions().nRegions()))
    {
        FatalErrorIn("crackerFvMesh::nCellsInRegion() const")
            << "region index is out of range"
            << abort(FatalError);
    }

    return (*nCellsInRegionPtr_)[regI];
}


const Foam::mapPolyMesh& Foam::crackerFvMesh::topoChangeMap() const
{
    return topoChangeMap_();
}


const Foam::vectorField& Foam::crackerFvMesh::globalCrackFaceCentres() const
{
    if (!globalCrackFaceCentresPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceCentresPtr_;
}


const Foam::scalarField& Foam::crackerFvMesh::globalCrackFaceSizes() const
{
    if (!globalCrackFaceSizesPtr_)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return *globalCrackFaceSizesPtr_;
}


const Foam::labelList& Foam::crackerFvMesh::globalCrackFaceAddressing() const
{
    if (!globalCrackFaceAddressingPtr_)
    {
        makeGlobalCrackFaceAddressing();
    }

    return *globalCrackFaceAddressingPtr_;
}


Foam::label Foam::crackerFvMesh::localCrackStart() const
{
    if (localCrackStart_ == -1)
    {
        makeGlobalCrackFaceCentresAndSizes();
    }

    return localCrackStart_;
}


Foam::label Foam::crackerFvMesh::globalCrackSize() const
{
    return globalCrackFaceCentres().size();
}

// ************************************************************************* //
