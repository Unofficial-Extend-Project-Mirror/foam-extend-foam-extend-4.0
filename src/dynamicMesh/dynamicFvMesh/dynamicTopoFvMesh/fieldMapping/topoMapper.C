/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    topoMapper

Description
    Implementation of topoMapper

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*----------------------------------------------------------------------------*/

#include "topoMapper.H"
#include "fluxCorrector.H"
#include "topoCellMapper.H"
#include "topoSurfaceMapper.H"
#include "topoBoundaryMeshMapper.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Store gradients prior to mesh reset
void topoMapper::storeGradients() const
{
    storeGradients<scalar>(sGrads_);
    storeGradients<vector>(vGrads_);

    if (fvMesh::debug)
    {
        Info<< "Registered volScalarFields: " << scalarGrads() << endl;
        Info<< "Registered volVectorFields: " << vectorGrads() << endl;
    }
}


//- Store geometric information
void topoMapper::storeGeometry() const
{
    // Wipe out existing information
    deleteDemandDrivenData(cellCentresPtr_);

    vectorField Cv(mesh_.cellCentres());
    vectorField Cf(mesh_.faceCentres());

    // Create and map the patch field values
    label nPatches = mesh_.boundary().size();

    // Create field parts
    PtrList<fvPatchField<vector> > volCentrePatches(nPatches);

    // Over-ride and set all patches to fixedValue
    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        volCentrePatches.set
        (
            patchI,
            new fixedValueFvPatchField<vector>
            (
                mesh_.boundary()[patchI],
                DimensionedField<vector, volMesh>::null()
            )
        );

        // Slice field to patch (forced assignment)
        volCentrePatches[patchI] ==
        (
            mesh_.boundaryMesh()[patchI].patchSlice(Cf)
        );
    }

    // Set the cell-centres pointer.
    cellCentresPtr_ =
    (
        new volVectorField
        (
            IOobject
            (
                "cellCentres",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimLength,
            SubField<vector>(Cv, mesh_.nCells()),
            volCentrePatches
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from mesh and dictionary
topoMapper::topoMapper
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    cellMap_(NULL),
    surfaceMap_(NULL),
    boundaryMap_(NULL),
    fluxCorrector_(fluxCorrector::New(mesh, dict)),
    cellCentresPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * *  //

topoMapper::~topoMapper()
{
    clear();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return reference to the mesh
const fvMesh& topoMapper::mesh() const
{
    return mesh_;
}


//- Return reference to objectRegistry storing fields.
const objectRegistry& topoMapper::thisDb() const
{
    return mesh_;
}


//- Set mapping information
void topoMapper::setMapper(const mapPolyMesh& mpm) const
{
    if
    (
        cellMap_.valid() ||
        surfaceMap_.valid() ||
        boundaryMap_.valid()
    )
    {
        FatalErrorIn
        (
            "void topoMapper::setMapper() const"
        ) << nl << " Mapper has already been set. "
          << abort(FatalError);
    }

    // Set pointers
    cellMap_.set(new topoCellMapper(mpm, *this));
    surfaceMap_.set(new topoSurfaceMapper(mpm, *this));
    boundaryMap_.set(new topoBoundaryMeshMapper(mesh(), mpm, *this));
}


//- Set face weighting information
void topoMapper::setFaceWeights
(
    const Xfer<List<scalarField> >& weights,
    const Xfer<List<vectorField> >& centres
) const
{
    faceWeights_.transfer(weights());
    faceCentres_.transfer(centres());
}


//- Set cell weighting information
void topoMapper::setCellWeights
(
    const Xfer<List<scalarField> >& weights,
    const Xfer<List<vectorField> >& centres
) const
{
    cellWeights_.transfer(weights());
    cellCentres_.transfer(centres());
}


//- Set cell / patch offset information
void topoMapper::setOffsets
(
    const labelList& cellSizes,
    const labelList& cellStarts,
    const labelList& faceSizes,
    const labelList& faceStarts,
    const labelListList& patchSizes,
    const labelListList& patchStarts
) const
{
    cellSizes_ = cellSizes;
    cellStarts_ = cellStarts;
    faceSizes_ = faceSizes;
    faceStarts_ = faceStarts;
    patchSizes_ = patchSizes;
    patchStarts_ = patchStarts;
}


//- Fetch face weights
const List<scalarField>& topoMapper::faceWeights() const
{
    return faceWeights_;
}


//- Fetch cell weights
const List<scalarField>& topoMapper::cellWeights() const
{
    return cellWeights_;
}


//- Fetch face centres
const List<vectorField>& topoMapper::faceCentres() const
{
    return faceCentres_;
}


//- Fetch cell centres
const List<vectorField>& topoMapper::cellCentres() const
{
    return cellCentres_;
}


//- Fetch cell sizes
const labelList& topoMapper::cellSizes() const
{
    return cellSizes_;
}


//- Fetch face sizes
const labelList& topoMapper::faceSizes() const
{
    return faceSizes_;
}


//- Fetch patch sizes
const labelListList& topoMapper::patchSizes() const
{
    return patchSizes_;
}


//- Fetch cell starts
const labelList& topoMapper::cellStarts() const
{
    return cellStarts_;
}


//- Fetch face starts
const labelList& topoMapper::faceStarts() const
{
    return faceStarts_;
}


//- Fetch patch starts
const labelListList& topoMapper::patchStarts() const
{
    return patchStarts_;
}


//- Store mesh information for the mapping stage
void topoMapper::storeMeshInformation() const
{
    // Store field-gradients
    storeGradients();

    // Store geometry
    storeGeometry();
}

//- Return non-const access to cell centres
volVectorField& topoMapper::volCentres() const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::volCentres() const"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellCentresPtr_;
}


//- Return stored cell centre information
const vectorField& topoMapper::internalCentres() const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::internalCentres() const"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellCentresPtr_;
}


//- Return stored patch centre information
const vectorField& topoMapper::patchCentres(const label i) const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::patchCentres"
            "(const label i) const"
        ) << nl << " Pointer has not been set. index: " << i
          << abort(FatalError);
    }

    return (*cellCentresPtr_).boundaryField()[i];
}


//- Return names of stored scalar gradients
const wordList topoMapper::scalarGrads() const
{
    return sGrads_.toc();
}


//- Return names of stored vector gradients
const wordList topoMapper::vectorGrads() const
{
    return vGrads_.toc();
}


//- Fetch the gradient field (template specialisation)
template <>
volVectorField& topoMapper::gradient(const word& name) const
{
    if (!sGrads_.found(name))
    {
        FatalErrorIn
        (
            "volVectorField& topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return sGrads_[name]();
}


//- Fetch the gradient field (template specialisation)
template <>
volTensorField& topoMapper::gradient(const word& name) const
{
    if (!vGrads_.found(name))
    {
        FatalErrorIn
        (
            "volTensorField& topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return vGrads_[name]();
}


//- Correct fluxes after topology changes, if required
void topoMapper::correctFluxes() const
{
    if (surfaceFluxCorrector().required())
    {
        // Supply a list of inserted faces for interpolation
        surfaceFluxCorrector().interpolateFluxes
        (
            surfaceMap().insertedObjectLabels()
        );

        // Update fluxes
        surfaceFluxCorrector().updateFluxes();
    }
}


//- Return volume mapper
const topoCellMapper& topoMapper::volMap() const
{
    if (!cellMap_.valid())
    {
        FatalErrorIn
        (
            "const topoCellMapper& topoMapper::volMap() const"
        ) << nl << " Volume mapper has not been set. "
          << abort(FatalError);
    }

    return cellMap_();
}


//- Return surface mapper
const topoSurfaceMapper& topoMapper::surfaceMap() const
{
    if (!surfaceMap_.valid())
    {
        FatalErrorIn
        (
            "const topoSurfaceMapper& topoMapper::surfaceMap() const"
        ) << nl << " Surface mapper has not been set. "
          << abort(FatalError);
    }

    return surfaceMap_();
}


//- Return boundary mapper
const topoBoundaryMeshMapper& topoMapper::boundaryMap() const
{
    if (!boundaryMap_.valid())
    {
        FatalErrorIn
        (
            "const topoBoundaryMeshMapper& topoMapper::boundaryMap() const"
        ) << nl << " Boundary mapper has not been set. "
          << abort(FatalError);
    }

    return boundaryMap_();
}


//- Return flux-corrector
const fluxCorrector& topoMapper::surfaceFluxCorrector() const
{
    if (!fluxCorrector_.valid())
    {
        FatalErrorIn
        (
            "const fluxCorrector& topoMapper::surfaceFluxCorrector() const"
        ) << nl << " fluxCorrector has not been set. "
          << abort(FatalError);
    }

    return fluxCorrector_();
}


//- Clear out member data
void topoMapper::clear() const
{
    // Clear out mappers
    cellMap_.clear();
    surfaceMap_.clear();
    boundaryMap_.clear();

    // Clear stored gradients
    sGrads_.clear();
    vGrads_.clear();

    // Wipe out geomtry information
    deleteDemandDrivenData(cellCentresPtr_);

    // Clear maps
    faceWeights_.clear();
    cellWeights_.clear();

    faceCentres_.clear();
    cellCentres_.clear();

    // Clear sizes / offsets
    cellSizes_.clear();
    cellStarts_.clear();

    faceSizes_.clear();
    faceStarts_.clear();

    patchSizes_.clear();
    patchStarts_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void topoMapper::operator=(const topoMapper& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "topoMapper::operator=(const topoMapper&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
