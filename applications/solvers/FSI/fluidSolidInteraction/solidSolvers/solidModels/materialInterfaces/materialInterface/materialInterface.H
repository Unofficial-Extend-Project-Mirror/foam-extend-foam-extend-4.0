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
    materialInterface

Description

SourceFiles
    materialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef MaterialInterface_H
#define MaterialInterface_H

#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "constitutiveModel.H"
#include "pointMesh.H"
#include "fvMatrices.H"
#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class materialInterface
:
    public regIOobject
{
    // Private data

        // Mesh
        const fvMesh& mesh_;

//         // Point mesh
//         const pointMesh pMesh_;

//         // Interface indicator point field
//         const pointScalarField materialInterface_;

    // Demand-driven data

        //- List of all faces at the interface
        mutable labelList* facesPtr_;

        //- Material sub-meshes (each sub-mesh contains only one material)
        mutable PtrList<fvMeshSubset> subMeshes_;

        //- Cell-to-point interpolation for each sub-mesh
        mutable PtrList<leastSquaresVolPointInterpolation> subMeshVolToPoint_;

        //- Number of materials for each point in the mesh
        mutable labelList* pointNumOfMaterialsPtr_;

        //- Isolated interface points
        mutable labelList* isolatedInterfacePointsPtr_;

    // Private Member Functions

        // Make demand-driven data

        //- Make global inter-faces addressing
        void makeFaces() const;

        //- Make material sub-meshes
        void makeSubMeshes() const;

        //- Make cell-to-point interpolators
        void makeSubMeshVolToPoint() const;

        //- Make materials number for points
        void makePointNumOfMaterials() const;

        //- Make isolated interface points
        void makeIsolatedInterfacePoints() const;

        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        materialInterface(const materialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const materialInterface&);

public:

    //- Runtime type information
    TypeName("materialInterface");

    // Constructors

        materialInterface(const fvMesh& mesh);

    // Destructor

        ~materialInterface();


    // Member Functions

        // Access

            //- Return reference to fvMesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return list of faces at the interface
            const labelList& faces() const;

            //- Return list of material sub-meshes
            const PtrList<fvMeshSubset>& subMeshes() const;

            //- Return list of cell-to-point interpolators
            const PtrList<leastSquaresVolPointInterpolation>&
            subMeshVolToPoint() const;

            //- Return material number for points
            const labelList& pointNumOfMaterials() const;

            //- Return isolated interface points
            const labelList& isolatedInterfacePoints() const;

            //- Name function is needed to disambiguate those inherited
            //  from regIOobject and dictionary
            const word& name() const
            {
                return regIOobject::name();
            }

        // Edit

            // Update interface after toplology change
            bool updateMesh(const mapPolyMesh& map)
            {
                clearOut();

                isolatedInterfacePoints();

                return true;
            }

            // Update interface after materials distribution change
            virtual bool update()
            {
                clearOut();

                isolatedInterfacePoints();

                return true;
            }

            //- Update displacement fields
            void updateDisplacement
            (
                const volVectorField& displacement,
                const vectorField& interfaceDisplacement,
                pointVectorField& pointDisplacement,
                PtrList<volVectorField>& subMeshDisplacement,
                PtrList<pointVectorField>& subMeshPointDisplacement
            );

            //- Calculate displacement gradients
            void updateDisplacementGradient
            (
                const volVectorField& displacement,
                const PtrList<volVectorField>& subMeshDisplacement,
                const PtrList<pointVectorField>& subMeshPointDisplacement,
                volTensorField& cellDisplacementGradient,
                surfaceTensorField& faceDisplacementGradient
            );

            //- Update displacement fields
            virtual void updateDisplacement
            (
                pointVectorField& pointDD
            )
            {}

            //- Calculate displacement gradients
            virtual void updateDisplacementGradient
            (
                volTensorField& gradDD,
                surfaceTensorField& gradDDf
            )
            {}

            //- Set properties to zero at the interface
            virtual void modifyProperty
            (
                surfaceScalarField& muf,
                const volScalarField& mu
            ) const;

            //- Interpolate properties from cells to facef
            //  Procedure at the interface depends on the particular
            //  material interface
            tmp<surfaceScalarField> interpolate
            (
                const volScalarField& mu
            ) const;

            //- Update total fields
            virtual void updateTotalFields()
            {}

            //- Correct discretization at the interface
            virtual void correct(fvVectorMatrix& DDEqn){}

            //- Correct traction at the interface
            virtual void correct(surfaceVectorField& t) const{}

            //- Calculate cell-centre gradient using Gauss method
            tmp<volTensorField> grad
            (
                const volVectorField& displacement,
                const vectorField& interfaceDisplacement
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
