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
    solidInterface

Description
    Special procedure to ensure accurate approximation of traction at discrete
    interfaces between cells of differing material properties.

    Tukovic et al. Int. J. Numer. Meth. Engng (2012) DOI: 10.1002/nme.4390

    This is the abstract base class and there are different specific
    implmentations depending on if there are thermal terms, nonlinear terms,
    non-orthogonal correction, orthotropic solver, etc.

Authors
    Zeljko Tukovic
    Aleksandar Karac
    Alojz Ivankovic
    Philip Cardff - reformatted to be allow runtime selection and added options

SourceFiles
    solidInterface.C
    newRheologyLaw.C

\*---------------------------------------------------------------------------*/

#ifndef solidInterface_H
#define solidInterface_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "constitutiveModel.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solidInterface Declaration
\*---------------------------------------------------------------------------*/

class solidInterface
{
    // Private data

        //- Name
        const word name_;

        const fvMesh& mesh_;

        const constitutiveModel& rheology_;


    // Demand-driven data

        //- Interface fvMeshSubMesh
        mutable fvMeshSubset* subMeshPtr_;

        //-
        mutable labelList* globalInterFacesPtr_;

        //-
        mutable labelList* localInterFacesPtr_;

        //-
        mutable vectorField* interfaceUPtr_;

        // depends on specific implementation
        //mutable volScalarField* muPtr_;
        //mutable volScalarField* lambdaPtr_;

        //- List of processor patches at the interface
        mutable labelList* processorPatchesPtr_;

        //- Processor faces at the interface
        mutable labelListList* processorPatchFacesPtr_;

        //- Interface displacement for processor patches
        mutable FieldField<Field, vector>* processorInterfaceUPtr_;

        //- Interface indicator
        mutable List<labelPair>* indicatorPtr_;

        //- whether to correct for non-orthogonality
        //Switch nonOrthogonalCorrection_;

        //- whether to include thermal stress terms
        //Switch thermalStress_;

        //- whether to include non-linear terms
        //Switch nonLinear_;

        // indicator field address
        // local address of interface face in globalInterFaces
        // -1 for non-interface faces
        mutable labelList* indicatorFieldMapPtr_;

        // Processor patches map from mesh to local processorPatches
        // -1 for non interface proc patches
        mutable labelList* processorPatchMapPtr_;

        // Processor faces map from mesh to local processorPatchFaces
        // local address into processorPatchFaces
        // -1 for non-interface faces
        mutable labelListList* processorPatchFacesMapPtr_;


    // Private Member Functions

        // Make demand-driven data

        //- Make sub-mesh
        void makeSubMesh() const;

        //- Make global inter-faces addressing
        void makeGlobalInterFaces() const;

        //- Make global inter-faces addressing
        void makeLocalInterFaces() const;

        //- Make interface displacement field
        void makeInterfaceDisplacement() const;

        //- Make processor patches
        void makeProcessorPatches() const;

        //- Make processor patch faces
        void makeProcessorPatchFaces() const;

        //- Make processor inter-faces displacement
        void makeProcessorInterfaceDisplacement() const;

        //- Make interface indicator
        void makeIndicator() const;

        //- Make interface indicator field
        //void makeIndicatorField() const;

        //- Make interface indicator field map
        void makeIndicatorFieldMap() const;

        //- Make processor interface indicator field map
        void makeProcessorPatchFacesMap() const;

        //- Clear all demand-driven data - philipc: make this public
        //void clearOut();

        //- Disallow copy construct
        solidInterface(const solidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const solidInterface&);


public:

    //- Runtime type information
    TypeName("solidInterface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidInterface,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const constitutiveModel& rheology
            ),
            (name, mesh, rheology)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<solidInterface> New
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Constructors

        //- Construct from dictionary
        solidInterface
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Destructor

      virtual ~solidInterface();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return const reference to the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const fvMeshSubset& subMesh() const;

        const labelList& globalInterFaces() const;

        const labelList& localInterFaces() const;

        vectorField& interfaceDisplacement();

        const vectorField& interfaceDisplacement() const;

        const labelList& processorPatches() const;

        const labelListList& processorPatchFaces() const;

        const FieldField<Field, vector>&
        processorInterfaceDisplacement() const;

        FieldField<Field, vector>& processorInterfaceDisplacement();

        virtual void correct(fvVectorMatrix& UEqn) = 0;

        virtual void modifyProperties
        (
            surfaceScalarField& s
        ) const;

        virtual void modifyProperties
        (
            surfaceSymmTensor4thOrderField& st
        ) const;

        virtual void modifyProperties
        (
            surfaceDiagTensorField& dt
        ) const;

        virtual void modifyProperties
        (
            surfaceScalarField& mu,
            surfaceScalarField& lambda
        ) const;

        virtual void modifyProperties
        (
            surfaceScalarField& mu,
            surfaceScalarField& lambda,
            surfaceScalarField& threeKalpha
        ) const;

        virtual void modifyProperties
        (
            surfaceSymmTensor4thOrderField& C,
            surfaceDiagTensorField& K
        ) const;

        const List<labelPair>& indicator() const;

        const labelList& indicatorFieldMap() const;

        const labelList& processorPatchMap() const;

        const labelListList& processorPatchFacesMap() const;

        virtual void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
