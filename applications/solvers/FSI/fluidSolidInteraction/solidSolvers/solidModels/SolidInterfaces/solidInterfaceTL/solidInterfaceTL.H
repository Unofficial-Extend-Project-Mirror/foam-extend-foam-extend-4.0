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
    solidInterfaceTL

Description

SourceFiles
    solidInterfaceTL.C

\*---------------------------------------------------------------------------*/

#ifndef SolidInterface_H
#define SolidInterface_H

#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "constitutiveModel.H"

#include "fvMatrices.H"
#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class solidInterfaceTL
:
    public regIOobject
{
    // Private data

        //- Displacement field
        const volVectorField& D_;

        //- Point displacement field
        const pointVectorField& pointD_;

//         //- Constitutive model
//         const constitutiveModel& rheology_;

    // Demand-driven data

        //- List of all faces at the interface
        mutable labelList* facesPtr_;

        //-
        mutable vectorField* displacementPtr_;

        //-
        mutable vectorField* tractionPtr_;

//         //-
//         mutable volScalarField* muPtr_;

//         //-
//         mutable volScalarField* lambdaPtr_;

        //- Material sub-meshes (each sub-mesh contains only one material)
        mutable PtrList<fvMeshSubset> subMeshes_;

        //- Cell-to-point interpolation for each sub-mesh
        mutable PtrList<leastSquaresVolPointInterpolation> subMeshVolToPoint_;

        //- Point displacement for each sub-mesh
        mutable PtrList<volVectorField> subMeshD_;

        //- Point displacement for each sub-mesh
        mutable PtrList<pointVectorField> subMeshPointD_;

        //- Number of materials for each point in the mesh
        mutable labelList* pointNumOfMaterialsPtr_;

    // Private Member Functions

        // Make demand-driven data

        //- Make global inter-faces addressing
        void makeFaces() const;

        //- Make interface displacement field
        void makeDisplacement() const;

        //- Make interface displacement field
        void makeTraction() const;

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Make cell-to-point interpolators
        void makeSubMeshVolToPoint() const;

        //- Make displacements field
        void makeSubMeshD() const;

        //- Make point displacements field
        void makeSubMeshPointD() const;

        //- Make materials number for points
        void makePointNumOfMaterials() const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        solidInterfaceTL(const solidInterfaceTL&);

        //- Disallow default bitwise assignment
        void operator=(const solidInterfaceTL&);

public:

    //- Runtime type information
    TypeName("solidInterfaceTL");

    // Constructors

        solidInterfaceTL
        (
            const volVectorField& D,
            const pointVectorField& pointD
//             const constitutiveModel& rheology
        );

    // Destructor

        virtual ~solidInterfaceTL();


    // Member Functions

        // Access

            //- Return list of faces at the interface
            const labelList& faces() const;

            //- Return interface displacement
            vectorField& displacement();

            //- Return interface displacement
            const vectorField& displacement() const;

            //- Return interface traction
            vectorField& traction();

            //- Return interface traction
            const vectorField& traction() const;

            //- Return list of material sub-meshes
            const PtrList<fvMeshSubset>& subMeshes() const;

            //- Return list of cell-to-point interpolators
            const PtrList<leastSquaresVolPointInterpolation>&
            subMeshVolToPoint() const;

            //- Return list of point displacement fields
            const PtrList<volVectorField>& subMeshD() const;

            //- Return list of point displacement fields
            PtrList<volVectorField>& subMeshD();

            //- Return list of point displacement fields
            const PtrList<pointVectorField>& subMeshPointD() const;

            //- Return list of point displacement fields
            PtrList<pointVectorField>& subMeshPointD();

            //- Return material number for points
            const labelList& pointNumOfMaterials() const;

            //- Name function is needed to disambiguate those inherited
            //  from regIOobject and dictionary
            const word& name() const
            {
                return regIOobject::name();
            }

        // Edit

            //- Correct discretization at the interface
            void correct(fvVectorMatrix& UEqn);

            //- Correct traction at the interface
            void correct(surfaceVectorField& t);

            //- Update displacement fields
            void updateDisplacement(pointVectorField& pointD);

            //- Calculate displacement gradients
            void updateDisplacementGradient
            (
                volTensorField& gradD,
                surfaceTensorField& gradDf
            );

            //- Set properties to zero at the interface
            void modifyProperties
            (
                surfaceScalarField& muf,
                surfaceScalarField& lambdaf
            ) const;

            //- Set properties to zero at the interface
            void modifyProperty
            (
                surfaceScalarField& muf
            ) const;

            //- Calculate cell-centre gradient using Gauss method
            tmp<volTensorField> grad(volVectorField& D) const;

            //-
            tmp<symmTensorField> sigmaA() const;

            //-
            tmp<symmTensorField> sigmaB() const;

            //- Interpolate from volField to pointField
            //  using least squares linear fit
            template<class Type>
            void volToPointInterpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                GeometricField<Type, pointPatchField, pointMesh>&
            ) const;

            //- Interpolate from volField to pointField
            //  using least squares linear fit
            template<class Type>
            void volToPointInterpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                const Field<Type>&,
                GeometricField<Type, pointPatchField, pointMesh>&
            ) const;


            //- ReadData function required for regIOobject read operation
            virtual bool readData(Istream&)
            {
                return false;
            }

            //- WriteData function required for regIOobject write operation
            virtual bool writeData(Ostream&) const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solidInterfaceTLVolPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
