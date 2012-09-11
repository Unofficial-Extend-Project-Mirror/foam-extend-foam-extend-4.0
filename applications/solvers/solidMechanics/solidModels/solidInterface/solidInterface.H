/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidInterface

Description
    special procedure to ensure accurate approximation of traction at discrete
    interfaces between cells of differing material properties

    Tukovic et al. Int. J. Numer. Meth. Engng (2012) DOI: 10.1002/nme.4390

SourceFiles
    solidInterface.C

Authors
    Zeljko Tukovic
    Aleksandar Karac
    Alojz Ivankovic

\*---------------------------------------------------------------------------*/

#ifndef SolidInterface_H
#define SolidInterface_H

#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "rheologyModel.H"

#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class solidInterface
{
    // Private data
    
        const fvMesh& mesh_;

        const rheologyModel& rheology_;

    // Demand-driven data

        //- Interface fvMeshSubMesh
        mutable fvMeshSubset* subMeshPtr_;

        //-
        mutable labelList* globalInterFacesPtr_;

        //-
        mutable labelList* localInterFacesPtr_;

        //-
        mutable vectorField* interfaceUPtr_;

        //-
        mutable volScalarField* muPtr_;

        //-
        mutable volScalarField* lambdaPtr_;

        //- List of processor patches at the interface
        mutable labelList* processorPatchesPtr_;

        //- Processor faces at the interface
        mutable labelListList* processorPatchFacesPtr_;

        //- Interface displacement for processor patches
        mutable FieldField<Field, vector>* processorInterfaceUPtr_;

        //- Interface indicator
        mutable List<labelPair>* indicatorPtr_;

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

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        solidInterface(const solidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const solidInterface&);

public:

    // Declare name of the class and it's debug switch
    ClassName("solidInterface");

    // Constructors
    solidInterface
    (
        const fvMesh& mesh,
        const rheologyModel& rheology
    );

    // Destructor
    ~solidInterface();  

    // Member Functions

    const fvMeshSubset& subMesh() const;    

    const labelList& globalInterFaces() const;

    const labelList& localInterFaces() const;

    vectorField& interfaceDisplacement();

    const vectorField& interfaceDisplacement() const;

    const labelList& processorPatches() const;

    const labelListList& processorPatchFaces() const;

    const FieldField<Field, vector>& processorInterfaceDisplacement() const;

    FieldField<Field, vector>& processorInterfaceDisplacement();

    void correct(fvVectorMatrix& UEqn);

    void modifyProperties
    (
        surfaceScalarField& muf, 
        surfaceScalarField& lambdaf
    ) const;

    tmp<volTensorField> grad(volVectorField& U) const;

    tmp<symmTensorField> sigmaA() const;

    tmp<symmTensorField> sigmaB() const;

    const List<labelPair>& indicator() const;

    void correctGrad(const volVectorField& U, volTensorField& gradU) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
