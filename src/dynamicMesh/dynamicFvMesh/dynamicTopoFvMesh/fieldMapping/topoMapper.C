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

#include "fvc.H"
#include "topoMapper.H"
#include "fluxCorrector.H"
#include "topoCellMapper.H"
#include "leastSquaresGrad.H"
#include "topoSurfaceMapper.H"
#include "topoBoundaryMeshMapper.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Store gradients of fields on the mesh prior to topology changes
template <class Type, class gradType>
void topoMapper::storeGradients
(
    HashTable<autoPtr<gradType> >& gradTable
) const
{
    // Define a few typedefs for convenience
    typedef typename outerProduct<vector, Type>::type gCmptType;
    typedef GeometricField<Type, fvPatchField, volMesh> volType;
    typedef GeometricField<gCmptType, fvPatchField, volMesh> gVolType;

    // Fetch all fields from registry
    HashTable<const volType*> fields
    (
        mesh_.objectRegistry::lookupClass<volType>()
    );

    forAllConstIter(typename HashTable<const volType*>, fields, fIter)
    {
        const volType& field = *fIter();

        // Compute the gradient.
        tmp<gVolType> tGrad;

        // If the fvSolution dictionary contains an entry,
        // use that, otherwise, default to leastSquares
        word gradName("grad(" + field.name() + ')');

        if (mesh_.schemesDict().subDict("gradSchemes").found(gradName))
        {
            tGrad = fvc::grad(field);
        }
        else
        {
            tGrad = fv::leastSquaresGrad<Type>(mesh_).grad(field);
        }

        // Make a new entry, but don't register the field.
        gradTable.insert
        (
            field.name(),
            autoPtr<gradType>
            (
                new gradType
                (
                    IOobject
                    (
                        tGrad().name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    tGrad()
                )
            )
        );
    }
}


//- Fetch the gradient field (template specialisation)
template <>
const volVectorField& topoMapper::gradient(const word& name) const
{
    if (!sGrads_.found(name))
    {
        FatalErrorIn
        (
            "const volVectorField& "
            "topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return sGrads_[name]();
}


//- Fetch the gradient field (template specialisation)
template <>
const volTensorField& topoMapper::gradient(const word& name) const
{
    if (!vGrads_.found(name))
    {
        FatalErrorIn
        (
            "const volTensorField& "
            "topoMapper::gradient(const word& name) const"
        ) << nl << " Gradient for: " << name
          << " has not been stored."
          << abort(FatalError);
    }

    return vGrads_[name]();
}


//- Store gradients prior to mesh reset
void topoMapper::storeGradients() const
{
    storeGradients<scalar>(sGrads_);
    storeGradients<vector>(vGrads_);

    if (fvMesh::debug)
    {
        Info << "Registered volScalarFields: " << endl;
        Info << sGrads_.toc() << endl;

        Info << "Registered volVectorFields: " << endl;
        Info << vGrads_.toc() << endl;
    }
}


//- Store geometric information
void topoMapper::storeGeometry() const
{
    if (cellCentresPtr_)
    {
        deleteDemandDrivenData(cellVolumesPtr_);
        deleteDemandDrivenData(cellCentresPtr_);

        patchAreasPtr_.clear();
        patchCentresPtr_.clear();
    }

    // Set the cell-volumes pointer.
    cellVolumesPtr_ = new scalarField(mesh_.cellVolumes());

    // Set the cell-centres pointer.
    cellCentresPtr_ = new vectorField(mesh_.cellCentres());

    // Set patch-areas
    patchAreasPtr_.setSize(mesh_.boundaryMesh().size());

    forAll(mesh_.boundaryMesh(), patchI)
    {
        patchAreasPtr_.set
        (
            patchI,
            new scalarField
            (
                mag(mesh_.boundaryMesh()[patchI].faceAreas())
            )
        );
    }

    // Set patch-centres.
    patchCentresPtr_.setSize(mesh_.boundaryMesh().size());

    forAll(mesh_.boundaryMesh(), patchI)
    {
        patchCentresPtr_.set
        (
            patchI,
            new vectorField
            (
                mesh_.boundaryMesh()[patchI].faceCentres()
            )
        );
    }
}


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


//- Store mesh information for the mapping stage
void topoMapper::storeMeshInformation() const
{
    // Store field-gradients
    storeGradients();

    // Store geometry
    storeGeometry();
}


//- Return stored cell volume information
const scalarField& topoMapper::cellVolumes() const
{
    if (!cellVolumesPtr_)
    {
        FatalErrorIn
        (
            "const scalarField& topoMapper::cellVolumes()"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellVolumesPtr_;
}


//- Return stored cell centre information
const vectorField& topoMapper::internalCentres() const
{
    if (!cellCentresPtr_)
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::internalCentres()"
        ) << nl << " Pointer has not been set. "
          << abort(FatalError);
    }

    return *cellCentresPtr_;
}


//- Return stored patch areas information
const scalarField& topoMapper::patchAreas(const label i) const
{
    if (!patchAreasPtr_.set(i))
    {
        FatalErrorIn
        (
            "const scalarField& topoMapper::patchAreas"
            "(const label i) const"
        ) << nl << " Pointer has not been set at index: " << i
          << abort(FatalError);
    }

    return patchAreasPtr_[i];
}


//- Return stored patch centre information
const vectorField& topoMapper::patchCentres(const label i) const
{
    if (!patchCentresPtr_.set(i))
    {
        FatalErrorIn
        (
            "const vectorField& topoMapper::patchCentres"
            "(const label i) const"
        ) << nl << " Pointer has not been set at index: " << i
          << abort(FatalError);
    }

    return patchCentresPtr_[i];
}


// Conservatively map all volFields in the registry
template <class Type>
void topoMapper::conservativeMapVolFields() const
{
    // Define a few typedefs for convenience
    typedef typename outerProduct<vector, Type>::type gCmptType;
    typedef GeometricField<Type, fvPatchField, volMesh> volType;
    typedef GeometricField<gCmptType, fvPatchField, volMesh> gradVolType;

    HashTable<const volType*> fields(mesh_.lookupClass<volType>());

    // Store old-times before mapping
    forAllIter(typename HashTable<const volType*>, fields, fIter)
    {
        volType& field = const_cast<volType&>(*fIter());

        field.storeOldTimes();
    }

    // Fetch internal/boundary mappers
    const topoCellMapper& fMap = volMap();
    const topoBoundaryMeshMapper& bMap = boundaryMap();

    // Now map all fields
    forAllIter(typename HashTable<const volType*>, fields, fIter)
    {
        volType& field = const_cast<volType&>(*fIter());

        if (fvMesh::debug)
        {
            Info << "Conservatively mapping "
                 << field.typeName
                 << ' ' << field.name()
                 << endl;
        }

        // Map the internal field
        fMap.mapInternalField
        (
            field.name(),
            gradient<gradVolType>(field.name()).internalField(),
            field.internalField()
        );

        // Map patch fields
        forAll(bMap, patchI)
        {
            bMap[patchI].mapPatchField
            (
                field.name(),
                field.boundaryField()[patchI]
            );
        }

        // Set the field instance
        field.instance() = field.mesh().thisDb().time().timeName();
    }
}


// Conservatively map all surfaceFields in the registry
template <class Type>
void topoMapper::conservativeMapSurfaceFields() const
{
    // Define a few typedefs for convenience
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfType;

    HashTable<const surfType*> fields(mesh_.lookupClass<surfType>());

    // Store old-times before mapping
    forAllIter(typename HashTable<const surfType*>, fields, fIter)
    {
        surfType& field = const_cast<surfType&>(*fIter());

        field.storeOldTimes();
    }

    // Fetch internal/boundary mappers
    const topoSurfaceMapper& fMap = surfaceMap();
    const topoBoundaryMeshMapper& bMap = boundaryMap();

    // Now map all fields
    forAllIter(typename HashTable<const surfType*>, fields, fIter)
    {
        surfType& field = const_cast<surfType&>(*fIter());

        if (fvMesh::debug)
        {
            Info << "Conservatively mapping "
                 << field.typeName
                 << ' ' << field.name()
                 << endl;
        }

        // Map the internal field
        fMap.mapInternalField
        (
            field.name(),
            field.internalField()
        );

        // Map patch fields
        forAll(bMap, patchI)
        {
            bMap[patchI].mapPatchField
            (
                field.name(),
                field.boundaryField()[patchI]
            );
        }

        // Set the field instance
        field.instance() = field.mesh().thisDb().time().timeName();
    }
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
            "const topoCellMapper& topoMapper::volMap()"
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
            "const topoSurfaceMapper& topoMapper::surfaceMap()"
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
            "const topoBoundaryMeshMapper& topoMapper::boundaryMap()"
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
            "const fluxCorrector& topoMapper::surfaceFluxCorrector()"
        ) << nl << " fluxCorrector has not been set. "
          << abort(FatalError);
    }

    return fluxCorrector_();
}


//- Clear out member data
void topoMapper::clear()
{
    // Clear out mappers
    cellMap_.clear();
    surfaceMap_.clear();
    boundaryMap_.clear();

    // Clear stored gradients
    sGrads_.clear();
    vGrads_.clear();

    // Wipe out geomtry information
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(cellCentresPtr_);
    patchAreasPtr_.clear();
    patchCentresPtr_.clear();

    // Clear maps
    faceWeights_.clear();
    cellWeights_.clear();

    faceCentres_.clear();
    cellCentres_.clear();
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
