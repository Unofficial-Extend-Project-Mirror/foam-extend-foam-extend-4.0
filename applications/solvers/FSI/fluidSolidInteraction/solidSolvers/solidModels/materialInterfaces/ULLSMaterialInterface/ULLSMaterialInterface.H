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
    ULLSMaterialInterface

Description

SourceFiles
    ULLSMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef ULLSMaterialInterface_H
#define ULLSMaterialInterface_H

#include "fvMesh.H"
// #include "fvMeshSubset.H"
#include "materialInterface.H"

#include "fvMatrices.H"
// #include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class ULLSMaterialInterface
:
    public materialInterface
{
    // Private data

        //- Displacement increment field
        const volVectorField& DD_;

        //- Point displacement increment field
        const pointVectorField& pointDD_;

        //- Displacement field
        const volVectorField& D_;

        //- Current time index
        label curTimeIndex_;

    // Demand-driven data

        //- Displacement increment at the interface
        mutable vectorField* displacementIncrementPtr_;

        //- Old displacement at the interface
        mutable vectorField* displacementPtr_;

//         //- Cauchy traction increment at the interface
//         mutable vectorField* tractionIncrementPtr_;

        //- Old Cauchy traction at the interface
        mutable vectorField* tractionPtr_;

        //- Relative deformation gradient at the interface
        mutable tensorField* ownRelFPtr_;

        //- Relative deformation gradient at the interface
        mutable tensorField* ngbRelFPtr_;

        //- Old deformation gradient at the interface
        mutable tensorField* ownFPtr_;

        //- Old deformation gradient at the interface
        mutable tensorField* ngbFPtr_;

        //- Old left Cauchy-Green deformation tensor at the interface
        mutable symmTensorField* ownBbarPtr_;

        //- Old left Cauchy-Green deformation tensor at the interface
        mutable symmTensorField* ngbBbarPtr_;


        //- Plastic strain tensor increment at the interface
        mutable symmTensorField* ownDEpsilonPPtr_;

        //- Plastic strain tensor increment at the interface
        mutable symmTensorField* ngbDEpsilonPPtr_;


        //- Plastic yield stress at the interface
        mutable scalarField* ownSigmaYPtr_;

        //- Plastic yield stress at the interface
        mutable scalarField* ngbSigmaYPtr_;


        //- Plastic yield stress increment at the interface
        mutable scalarField* ownDSigmaYPtr_;

        //- Plastic yield stress increment at the interface
        mutable scalarField* ngbDSigmaYPtr_;


        //- Plastic yield stress at the interface
        mutable scalarField* ownEpsilonPEqPtr_;

        //- Plastic yield stress at the interface
        mutable scalarField* ngbEpsilonPEqPtr_;


        //- Plastic yield stress at the interface
        mutable scalarField* ownDEpsilonPEqPtr_;

        //- Plastic yield stress at the interface
        mutable scalarField* ngbDEpsilonPEqPtr_;


        //- Point displacement increment for each sub-mesh
        mutable PtrList<volVectorField> subMeshDD_;

        //- Point displacement increment for each sub-mesh
        mutable PtrList<pointVectorField> subMeshPointDD_;

    // Private Member Functions

        // Make demand-driven data

        //- Make interface displacement field
        void makeDisplacementIncrement() const;

        //- Make interface displacement field
        void makeDisplacement() const;

        //- Make interface displacement field
        void makeTractionIncrement() const;

        //- Make interface displacement field
        void makeTraction() const;

        //- Make relative deformation gradient at the interface
        void makeRelF() const;

        //- Make old deformation gradient at the interface
        void makeF() const;

        //- Make old left Cauchy-Green deformation tensor at the interface
        void makeBbar() const;

        //- Make plastic strain increment at the interface
        void makeDEpsilonP() const;

        //- Make yield stress at the interface
        void makeSigmaY() const;

        //- Make yield stress increment at the interface
        void makeDSigmaY() const;

        //- Make equivalent plastic strain at the interface
        void makeEpsilonPEq() const;

        //- Make icrement of equivalent plastic strain at the interface
        void makeDEpsilonPEq() const;

        //- Make displacements field
        void makeSubMeshDD() const;

        //- Make point displacements field
        void makeSubMeshPointDD() const;

        //- Cauchy traction at the owner and neighbour side of the interface
        tmp<vectorField> cauchyTraction
        (
            const vector& m,
            const label faceI,
            const volTensorField& gradDD,
            const surfaceTensorField& gradDDf,
            const volScalarField& mu,
            const volScalarField& lambda,
            const bool plasticityActive = false
        ) const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        ULLSMaterialInterface(const ULLSMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const ULLSMaterialInterface&);

public:

    //- Runtime type information
    TypeName("ULLSMaterialInterface");

    // Constructors

        ULLSMaterialInterface
        (
            const volVectorField& DD,
            const pointVectorField& pointDD,
            const volVectorField& D
        );

    // Destructor

        virtual ~ULLSMaterialInterface();


    // Member Functions

        // Access

            //- Return reference to displacement increment
            const volVectorField& DD() const
            {
                return DD_;
            }

            //- Return interface displacement increment
            vectorField& displacementIncrement();

            //- Return interface displacement increment
            const vectorField& displacementIncrement() const;

            //- Return interface displacement
            vectorField& displacement();

            //- Return interface displacement
            const vectorField& displacement() const;

//             //- Return interface traction increment
//             vectorField& tractionIncrement();

//             //- Return interface traction increment
//             const vectorField& tractionIncrement() const;

            //- Return interface traction
            vectorField& traction();

            //- Return interface traction
            const vectorField& traction() const;

            //- Return relative deformation gradient at the interface
            tensorField& ownRelF();

            //- Return relative deformation gradient at the interface
            const tensorField& ownRelF() const;

            //- Return relative deformation gradient at the interface
            tensorField& ngbRelF();

            //- Return relative deformation gradient at the interface
            const tensorField& ngbRelF() const;

            //- Return old deformation gradient at the interface
            tensorField& ownF();

            //- Return old deformation gradient at the interface
            const tensorField& ownF() const;

            //- Return old deformation gradient at the interface
            tensorField& ngbF();

            //- Return old deformation gradient at the interface
            const tensorField& ngbF() const;

            //- Return old left Cauchy-Green deformation tensor
            symmTensorField& ownBbar();

            //- Return old left Cauchy-Green deformation tensor
            const symmTensorField& ownBbar() const;

            //- Return old left Cauchy-Green deformation tensor
            symmTensorField& ngbBbar();

            //- Return old left Cauchy-Green deformation tensor
            const symmTensorField& ngbBbar() const;


            //- Return plastic strain increment at the interface
            symmTensorField& ownDEpsilonP();

            //- Return plastic strain increment at the interface
            const symmTensorField& ownDEpsilonP() const;

            //- Return plastic strain increment at the interface
            symmTensorField& ngbDEpsilonP();

            //- Return plastic strain increment at the interface
            const symmTensorField& ngbDEpsilonP() const;


            //- Return yeald stress at the interface
            scalarField& ownSigmaY();

            //- Return yeald stress at the interface
            const scalarField& ownSigmaY() const ;

            //- Return yeald stress at the interface
            scalarField& ngbSigmaY();

            //- Return yeald stress at the interface
            const scalarField& ngbSigmaY() const ;


            //- Return yeald stress increment at the interface
            scalarField& ownDSigmaY();

            //- Return yeald stress increment at the interface
            const scalarField& ownDSigmaY() const ;

            //- Return yeald stress increment at the interface
            scalarField& ngbDSigmaY();

            //- Return yeald stress increment at the interface
            const scalarField& ngbDSigmaY() const ;


            //- Return ... at the interface
            scalarField& ownEpsilonPEq();

            //- Return ... at the interface
            const scalarField& ownEpsilonPEq() const ;

            //- Return ... at the interface
            scalarField& ngbEpsilonPEq();

            //- Return ... at the interface
            const scalarField& ngbEpsilonPEq() const ;


            //- Return ... at the interface
            scalarField& ownDEpsilonPEq();

            //- Return ... at the interface
            const scalarField& ownDEpsilonPEq() const ;

            //- Return ... at the interface
            scalarField& ngbDEpsilonPEq();

            //- Return ... at the interface
            const scalarField& ngbDEpsilonPEq() const ;


            //- Return list of point displacement fields
            const PtrList<volVectorField>& subMeshDD() const;

            //- Return list of point displacement fields
            PtrList<volVectorField>& subMeshDD();

            //- Return list of point displacement fields
            const PtrList<pointVectorField>& subMeshPointDD() const;

            //- Return list of point displacement fields
            PtrList<pointVectorField>& subMeshPointDD();

        // Edit

            //- Update interface after toplogoy change
            bool updateMesh(const mapPolyMesh&);

            //- Update interface after materials distribution change
            virtual bool update();

            //- Correct discretization at the interface
            virtual void correct(fvVectorMatrix& DDEqn);

            //- Correct traction at the interface
            virtual void correct(surfaceVectorField& t) const;

            //- Update interface displacement increment field
            void updateInterfaceDisplacementIncrement();

            //-  Update normal component of interface displacement gradient
            void correctInterfaceGradient(surfaceTensorField& gradDDf);

            //- Update displacement fields
            void updateDisplacementIncrement(pointVectorField& pointDD);

            //- Calculate displacement gradients
            void updateDisplacementIncrementGradient
            (
                volTensorField& gradDD,
                surfaceTensorField& gradDDf
            );

            //- Update displacement fields
            virtual void updateDisplacement
            (
                pointVectorField& pointDD
            )
            {
                updateDisplacementIncrement(pointDD);
            }

            //- Calculate displacement gradients
            virtual void updateDisplacementGradient
            (
                volTensorField& gradDD,
                surfaceTensorField& gradDDf
            )
            {
                updateDisplacementIncrementGradient(gradDD, gradDDf);
            }

            //- Update total fields
            virtual void updateTotalFields();

            //- Set properties to zero at the interface
            virtual void modifyProperty
            (
                surfaceScalarField& muf,
                const volScalarField& mu
            ) const;

            //- Return trial left Cauchy-Green deformation tensor
            //  at the owner and neighbour side of the interface
            tmp<symmTensorField> trialBbar
            (
                const label faceI,
                const volTensorField& gradDD,
                const surfaceTensorField& gradDDf
            ) const;

            //- Return determinant of deformation gradient
            //  at the owner and neighbour side of the interface
            tmp<scalarField> J
            (
                const label faceI,
                const volTensorField& gradDD,
                const surfaceTensorField& gradDDf
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
