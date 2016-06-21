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
    solidInterfaceITL

Description

SourceFiles
    solidInterfaceITL.C

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

class solidInterfaceITL
:
    public regIOobject
{
    // Private data

        //- Displacement increment field
        const volVectorField& DD_;

        //- Displacement field
        const volVectorField& D_;

        //- Current time index
        label curTimeIndex_;

    // Demand-driven data

        //- List of all faces at the interface
        mutable labelList* facesPtr_;

        //- Displacement increment at the interface
        mutable vectorField* displacementIncrementPtr_;

        //- Old displacement at the interface
        mutable vectorField* displacementPtr_;

        //- Cauchy traction increment at the interface
        mutable vectorField* tractionIncrementPtr_;

        //- Old Cauchy traction at the interface
        mutable vectorField* tractionPtr_;

        //- Second Piola-Kirchhoff stress increment at the interface
        mutable symmTensorField* ownDSigmafPtr_;

        //- Second Piola-Kirchhoff stress increment at the interface
        mutable symmTensorField* ngbDSigmafPtr_;

        //- Old second Piola-Kirchhoff stress at the interface
        mutable symmTensorField* ownSigmafPtr_;

        //- Old second Piola-Kirchhoff stress at the interface
        mutable symmTensorField* ngbSigmafPtr_;

        //- Material sub-meshes (each sub-mesh contains only one material)
        mutable PtrList<fvMeshSubset> subMeshes_;

        //- Cell-to-point interpolation for each sub-mesh
        mutable PtrList<leastSquaresVolPointInterpolation> subMeshVolToPoint_;

        //- Number of materials for each point in the mesh
        mutable labelList* pointNumOfMaterialsPtr_;

    // Private Member Functions

        // Make demand-driven data

        //- Make global inter-faces addressing
        void makeFaces() const;

        //- Make interface displacement field
        void makeDisplacementIncrement() const;

        //- Make interface displacement field
        void makeDisplacement() const;

        //- Make interface displacement field
        void makeTractionIncrement() const;

        //- Make interface displacement field
        void makeTraction() const;

        //- Make second Piola-Kirchhoff stress increment at the interface
        void makeDSigmaf() const;

        //- Make old second Piola-Kirchhoff stress at the interface
        void makeSigmaf() const;

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Make cell-to-point interpolators
        void makeSubMeshVolToPoint() const;

        //- Make materials number for points
        void makePointNumOfMaterials() const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        solidInterfaceITL(const solidInterfaceITL&);

        //- Disallow default bitwise assignment
        void operator=(const solidInterfaceITL&);

public:

    //- Runtime type information
    TypeName("solidInterfaceITL");

    // Constructors

        solidInterfaceITL
        (
            const volVectorField& DD,
            const volVectorField& D
        );

    // Destructor

        virtual ~solidInterfaceITL();


    // Member Functions

        // Access

            //- Return list of faces at the interface
            const labelList& faces() const;

            //- Return interface displacement increment
            vectorField& displacementIncrement();

            //- Return interface displacement increment
            const vectorField& displacementIncrement() const;

            //- Return interface displacement
            vectorField& displacement();

            //- Return interface displacement
            const vectorField& displacement() const;

            //- Return interface traction increment
            vectorField& tractionIncrement();

            //- Return interface traction increment
            const vectorField& tractionIncrement() const;

            //- Return interface traction
            vectorField& traction();

            //- Return interface traction
            const vectorField& traction() const;

            //- Return old second Piola-Kirchhoff stress
            symmTensorField& ownSigmaf();

            //- Return old second Piola-Kirchhoff stress
            const symmTensorField& ownSigmaf() const;

            //- Return old second Piola-Kirchhoff stress
            symmTensorField& ngbSigmaf();

            //- Return old second Piola-Kirchhoff stress
            const symmTensorField& ngbSigmaf() const;

            //- Return old second Piola-Kirchhoff stress
            symmTensorField& ownDSigmaf();

            //- Return old second Piola-Kirchhoff stress
            const symmTensorField& ownDSigmaf() const;

            //- Return old second Piola-Kirchhoff stress
            symmTensorField& ngbDSigmaf();

            //- Return old second Piola-Kirchhoff stress
            const symmTensorField& ngbDSigmaf() const;

            //- Return list of material sub-meshes
            const PtrList<fvMeshSubset>& subMeshes() const;

            //- Return list of cell-to-point interpolators
            const PtrList<leastSquaresVolPointInterpolation>&
            subMeshVolToPoint() const;

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
            void correct(fvVectorMatrix& DDEqn);

            //- Correct traction at the interface
            void correct(surfaceVectorField& t);

            //- Set properties to zero at the interface
            void modifyProperties
            (
                surfaceScalarField& muf,
                surfaceScalarField& lambdaf
            ) const;

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
#   include "solidInterfaceITLVolPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
