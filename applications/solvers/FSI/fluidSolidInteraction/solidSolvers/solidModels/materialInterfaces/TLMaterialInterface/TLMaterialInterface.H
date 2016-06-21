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
    TLMaterialInterface

Description

SourceFiles
    TLMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef TLMaterialInterface_H
#define TLMaterialInterface_H

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

class TLMaterialInterface
:
    public materialInterface
{
    // Private data

        //- Displacement field
        const volVectorField& D_;

        //- Point displacement field
        const pointVectorField& pointD_;

    // Demand-driven data

        //- Displacement at the material interface faces
        mutable vectorField* displacementPtr_;

        //- Traction a the material interface faces
        mutable vectorField* tractionPtr_;

        //- Point displacement for each sub-mesh
        mutable PtrList<volVectorField> subMeshD_;

        //- Point displacement for each sub-mesh
        mutable PtrList<pointVectorField> subMeshPointD_;

    // Private Member Functions

        // Make demand-driven data

        //- Make interface displacement field
        void makeDisplacement() const;

        //- Make interface displacement field
        void makeTraction() const;

        //- Make displacements field
        void makeSubMeshD() const;

        //- Make point displacements field
        void makeSubMeshPointD() const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        TLMaterialInterface(const TLMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const TLMaterialInterface&);

public:

    //- Runtime type information
    TypeName("TLMaterialInterface");

    // Constructors

        TLMaterialInterface
        (
            const volVectorField& D,
            const pointVectorField& pointD
        );

    // Destructor

        ~TLMaterialInterface();


    // Member Functions

        // Access

            //- Return interface displacement
            vectorField& displacement();

            //- Return interface displacement
            const vectorField& displacement() const;

            //- Return interface traction
            vectorField& traction();

            //- Return interface traction
            const vectorField& traction() const;

            //- Return list of point displacement fields
            const PtrList<volVectorField>& subMeshD() const;

            //- Return list of point displacement fields
            PtrList<volVectorField>& subMeshD();

            //- Return list of point displacement fields
            const PtrList<pointVectorField>& subMeshPointD() const;

            //- Return list of point displacement fields
            PtrList<pointVectorField>& subMeshPointD();

        // Edit

            //- Update interface after materials distribution change
            virtual bool update();

            //- Correct discretization at the interface
            virtual void correct(fvVectorMatrix& UEqn);

            //- Correct traction at the interface
            virtual void correct(surfaceVectorField& t) const;

            //- Update displacement fields
            virtual void updateDisplacement(pointVectorField& pointD);

            //- Calculate displacement gradients
            virtual void updateDisplacementGradient
            (
                volTensorField& gradD,
                surfaceTensorField& gradDf
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
