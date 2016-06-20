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

#include "Tuple2.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "IOmanip.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fluxCorrector;
class topoCellMapper;
class topoSurfaceMapper;
class topoBoundaryMeshMapper;

/*---------------------------------------------------------------------------*\
                         Class topoMapper Declaration
\*---------------------------------------------------------------------------*/

class topoMapper
{
    // Private typedefs

        typedef Tuple2<word, label> GradientMap;
        typedef HashTable<GradientMap> GradientTable;

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

        // Index map for gradients
        mutable GradientTable gradTable_;

        // Stored gradients for mapping
        mutable PtrList<volVectorField> sGradPtrs_;
        mutable PtrList<volTensorField> vGradPtrs_;

        //- Geometric information on the old mesh
        mutable scalarField* cellVolumesPtr_;
        mutable volVectorField* cellCentresPtr_;

        //- Intersection weights
        mutable List<scalarField> faceWeights_;
        mutable List<scalarField> cellWeights_;

        //- Intersection centres
        mutable List<vectorField> faceCentres_;
        mutable List<vectorField> cellCentres_;

        //- Sizes / starts for mapping
        mutable labelList cellSizes_;
        mutable labelList cellStarts_;
        mutable labelList faceSizes_;
        mutable labelList faceStarts_;
        mutable labelListList patchSizes_;
        mutable labelListList patchStarts_;

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
            GradientTable& gradTable,
            PtrList<gradType>& gradList
        ) const;

        //- Store gradients prior to mesh reset
        void storeGradients() const;

        //- Set geometric information
        void storeGeometry() const;

public:

    // Constructors

        //- Construct from mesh and dictionary
        topoMapper(const fvMesh& mesh, const dictionary& dict);

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

        //- Set cell / patch offset information
        void setOffsets
        (
            const labelList& cellSizes,
            const labelList& cellStarts,
            const labelList& faceSizes,
            const labelList& faceStarts,
            const labelListList& patchSizes,
            const labelListList& patchStarts
        ) const;

        //- Fetch face weights
        const List<scalarField>& faceWeights() const;

        //- Fetch cell weights
        const List<scalarField>& cellWeights() const;

        //- Fetch face centres
        const List<vectorField>& faceCentres() const;

        //- Fetch cell centres
        const List<vectorField>& cellCentres() const;

        //- Fetch cell sizes
        const labelList& cellSizes() const;

        //- Fetch face sizes
        const labelList& faceSizes() const;

        //- Fetch patch sizes
        const labelListList& patchSizes() const;

        //- Fetch cell starts
        const labelList& cellStarts() const;

        //- Fetch face starts
        const labelList& faceStarts() const;

        //- Fetch patch starts
        const labelListList& patchStarts() const;

        //- Store mesh information for the mapping stage
        void storeMeshInformation() const;

        //- Deregister gradient fields and centres,
        //  but retain for mapping
        void deregisterMeshInformation() const;

        //- Return non-const access to cell centres
        volVectorField& volCentres() const;

        //- Return non-const access to cell volumes
        scalarField& internalVolumes() const;

        //- Return stored cell centre information
        const vectorField& internalCentres() const;

        //- Return stored patch centre information
        const vectorField& patchCentres(const label i) const;

        //- Return names of stored gradients
        const wordList gradientTable() const;

        //- Fetch the gradient field
        template <class Type>
        Type& gradient(const word& name) const;

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
        void clear() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "topoMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
