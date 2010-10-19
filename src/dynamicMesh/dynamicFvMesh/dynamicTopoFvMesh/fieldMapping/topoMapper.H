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
    Class holds all necessary information for mapping fields associated with
    dynamicTopoFvMesh and fvMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoMapper_H
#define topoMapper_H

#include "fluxCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class topoCellMapper;
class topoSurfaceMapper;
class topoBoundaryMeshMapper;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class topoMapper Declaration
\*---------------------------------------------------------------------------*/

class topoMapper
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Reference to the options dictionary
        const dictionary& dict_;

    // Demand-driven private data

        //- Cell mapper
        mutable autoPtr<topoCellMapper> cellMap_;

        //- Surface mapper
        mutable autoPtr<topoSurfaceMapper> surfaceMap_;

        //- Boundary mapper
        mutable autoPtr<topoBoundaryMeshMapper> boundaryMap_;

        //- Flux corrector
        mutable autoPtr<fluxCorrector> fluxCorrector_;

        // Stored gradients for mapping
        mutable HashTable<autoPtr<volVectorField> > sGrads_;
        mutable HashTable<autoPtr<volTensorField> > vGrads_;

        //- Geometric information on the old mesh
        mutable scalarField* cellVolumesPtr_;
        mutable vectorField* cellCentresPtr_;
        mutable PtrList<scalarField> patchAreasPtr_;
        mutable PtrList<vectorField> patchCentresPtr_;

        //- Intersection volume weights
        mutable List<scalarField> faceWeights_;
        mutable List<scalarField> cellWeights_;

        //- Intersection centre weights
        mutable List<vectorField> faceCentres_;
        mutable List<vectorField> cellCentres_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoMapper(const topoMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoMapper&);

        // Store gradients of volFields on the mesh
        // prior to topology changes
        template <class Type, class gradType>
        void storeGradients
        (
            HashTable<autoPtr<gradType> >& gradTable
        ) const;

        //- Fetch the gradient field
        template <class Type>
        const Type& gradient(const word& name) const;

        //- Store gradients prior to mesh reset
        void storeGradients() const;

        //- Set geometric information
        void storeGeometry() const;

public:

    // Constructors

        //- Construct from mesh and dictionary
        topoMapper(const fvMesh& mesh, const dictionary& dict)
        :
            mesh_(mesh),
            dict_(dict),
            cellMap_(NULL),
            surfaceMap_(NULL),
            boundaryMap_(NULL),
            fluxCorrector_(fluxCorrector::New(mesh, dict)),
            cellVolumesPtr_(NULL),
            cellCentresPtr_(NULL),
            patchAreasPtr_(mesh.boundary().size()),
            patchCentresPtr_(mesh.boundary().size())
        {}

    // Destructor

        ~topoMapper();

    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const;

        //- Return reference to objectRegistry storing fields.
        const objectRegistry& thisDb() const;

        //- Set mapping information
        void setMapper(const mapPolyMesh& mpm) const;

        //- Set face weighting information
        void setFaceWeights
        (
            const Xfer<List<scalarField> >& weights,
            const Xfer<List<vectorField> >& centres
        ) const;

        //- Set cell weighting information
        void setCellWeights
        (
            const Xfer<List<scalarField> >& weights,
            const Xfer<List<vectorField> >& centres
        ) const;

        //- Fetch face weights
        const List<scalarField>& faceWeights() const;

        //- Fetch cell weights
        const List<scalarField>& cellWeights() const;

        //- Fetch face centres
        const List<vectorField>& faceCentres() const;

        //- Fetch cell centres
        const List<vectorField>& cellCentres() const;

        //- Store mesh information for the mapping stage
        void storeMeshInformation() const;

        //- Return stored cell volume information
        const scalarField& cellVolumes() const;

        //- Return stored cell centre information
        const vectorField& internalCentres() const;

        //- Return stored patch area information
        const scalarField& patchAreas(const label i) const;

        //- Return stored patch centre information
        const vectorField& patchCentres(const label i) const;

        // Conservatively map all volFields in the registry
        template <class Type>
        void conservativeMapVolFields() const;

        // Conservatively map all surfaceFields in the registry
        template <class Type>
        void conservativeMapSurfaceFields() const;

        //- Correct fluxes after topology change
        void correctFluxes() const;

        //- Return volume mapper
        const topoCellMapper& volMap() const;

        //- Return surface mapper
        const topoSurfaceMapper& surfaceMap() const;

        //- Return boundary mapper
        const topoBoundaryMeshMapper& boundaryMap() const;

        //- Return flux-corrector
        const fluxCorrector& surfaceFluxCorrector() const;

        //- Clear out member data
        void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "topoMapper.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
