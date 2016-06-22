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
    Foam::leastSquaresVolPointInterpolation

Description
    Foam::leastSquaresVolPointInterpolation

SourceFiles
    leastSquaresVolPointInterpolation.C
    leastSquaresVolPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVolPointInterpolation_H
#define leastSquaresVolPointInterpolation_H

#include "MeshObject.H"
#include "PtrList.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "className.H"
#include "scalarMatrices.H"
#include "FieldField.H"
#include "Map.H"
#include "polyPatch.H"
#include "Tuple2.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class leastSquaresVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresVolPointInterpolation
:
    public MeshObject<fvMesh, leastSquaresVolPointInterpolation>
{
    // Private data

        //- Point boundary face addressig
        mutable labelListList* pointBndFacesPtr_;

        //- Point cyclic boundary face addressig
        mutable labelListList* pointCyclicFacesPtr_;

        //- Point boundary face addressig
        mutable labelListList* pointCyclicBndFacesPtr_;

        //- Point cyclic ggi boundary face addressig
        mutable labelListList* pointCyclicGgiFacesPtr_;

        //- Point cyclic ggi boundary face addressig
        mutable labelListList* pointCyclicGgiBndFacesPtr_;

        //- Point boundary face addressig
        mutable labelListList* pointProcFacesPtr_;

        //- Axis edges for axisymmetric cases
        mutable labelList* axisEdgesPtr_;

        //- Point ngb processor boundary face centres
        mutable Map<labelList>* pointAxisEdgesPtr_;

//         //- Point ngb processor boundary face centres
//         mutable Map<Field<vector> >* pointNgbProcBndFaceCentresPtr_;

        //- Global point ngb processor boundary face centres
        mutable Map<Field<vector> >* globalPointNgbProcBndFaceCentresPtr_;

        //- Global point ngb processor cell centres
        mutable Map<Field<vector> >* globalPointNgbProcCellCentresPtr_;

        //- Locel boundary faces needed by neighbour processors
        mutable labelListList* procBndFacesPtr_;

        //- Centres of boundary faces from neigbour processores
        mutable FieldField<Field, vector>* procBndFaceCentresPtr_;

//         //- Point ngb processor boundary faces addressing
//         mutable Map<List<labelPair> >* pointProcBndFacesPtr_;

        //- Point ngb processor boundary faces addressing
        mutable List<List<labelPair> >* pointProcBndFacesPtr_;

        //- Locel processor cells needed by neighbour processors
        mutable labelListList* procCellsPtr_;

        //- Point ngb processor cell addressing
        mutable Map<List<labelPair> >* pointProcCellsPtr_;

        //- Centres of cells from neigbour processores
        mutable FieldField<Field, vector>* procCellCentresPtr_;

        //- Weighting factors
        mutable FieldField<Field, scalar>* weightsPtr_;

        //- Local origin
        mutable Field<vector>* originsPtr_;

        //- Mirror plane normals and transformation tensors
        mutable List<Tuple2<vector, tensor> >* mirrorPlaneTransformationPtr_;
//         mutable Map<Tuple2<vector, tensor> >* mirrorPlaneTransformationPtr_;

        //- Interpolation matrices
        mutable PtrList<scalarRectangularMatrix> invLsMatrices_;

        //- Interpolatin stencil size
        mutable scalarField* refLPtr_;

    // Private member functions

        //- Make point boundary face addressing
        void makePointFaces() const;

        //- Make axis edges
        void makeAxisEdges() const;

        //- Make point axis edges addressing
        void makePointAxisEdges() const;

//         //- Make point ngb processor boundary face centres
//         void makePointNgbProcBndFaceCentres() const;

        //- Make global point ngb processor cell centres
        void makeGlobalPointNgbProcBndFaceCentres() const;

        //- Make global point ngb processor cell centres
        void makeGlobalPointNgbProcCellCentres() const;

        //- Make local boundary faces needed by neighbour processors
        void makeProcBndFaces() const;

        //- Make boundary face centres from neighbour processors
        void makeProcBndFaceCentres() const;

        //- Make cells needed by neighbour processors
        void makeProcCells() const;

        //- Make cell centres from neighbour processors
        void makeProcCellCentres() const;

        //- Make weights
        void makeWeights() const;

        //- Make local origins
        void makeOrigins() const;

        //- Make interpolation matrix
        void makeInvLsMatrices() const;

        //- Make mirror plane noramals and transformation tensors
        void makeMirrorPlaneTransformation() const;

        //- Disallow default bitwise copy construct
        leastSquaresVolPointInterpolation
        (
            const leastSquaresVolPointInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresVolPointInterpolation&);

public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresVolPointInterpolation");


    // Constructors

        //- Constructor given fvMesh.  pointMesh will be created or
        //  looked up from objectRegistry
        explicit leastSquaresVolPointInterpolation(const fvMesh&);


    // Destructor

        virtual ~leastSquaresVolPointInterpolation();


    // Member functions

        // Access

            //- Get point boundary face addressing
            const labelListList& pointBndFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicBndFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicGgiFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicGgiBndFaces() const;

            //- Get point processor boundary face addressing
            const labelListList& pointProcFaces() const;

            //- Get axis edges
            const labelList& axisEdges() const;

            //- Get point axis edges addressing
            const Map<labelList>& pointAxisEdges() const;

//             //- Get point processor boundary face addressing
//             const Map<Field<vector> >& pointNgbProcBndFaceCentres() const;

            //- Get global point ngb processor boundary face centres
            const Map<Field<vector> >&
            globalPointNgbProcBndFaceCentres() const;

            //- Get global point ngb processor cell centres
            const Map<Field<vector> >& globalPointNgbProcCellCentres() const;

            //- Get list of boundary faces needed by neighbour processores
            const labelListList& procBndFaces() const;

            //- Get centres of boundary faces from neighbour processors
            const FieldField<Field, vector>& procBndFaceCentres() const;

            //- Get neighbor processor cell point addresing
            const List<List<labelPair> >& pointProcBndFaces() const;
//             const Map<List<labelPair> >& pointProcBndFaces() const;

            //- Get list of cells needed by neighbour processores
            const labelListList& procCells() const;

            //- Get neighbor processor cell point addresing
            const Map<List<labelPair> >& pointProcCells() const;

            //- Get centres of cells from neighbour processors
            const FieldField<Field, vector>& procCellCentres() const;

            //- Get local origins
            const FieldField<Field, scalar>& weights() const;

            //- Get local origins
            const Field<vector>& origins() const;

            //- Get interpolation matrices
            const PtrList<scalarRectangularMatrix>& invLsMatrices() const;

            //- Get interpolation stencil size
            const scalarField& refL() const;

            //- Get mirror plane normals
            const List<Tuple2<vector, tensor> >&
            mirrorPlaneTransformation() const;
//             const Map<Tuple2<vector, tensor> >&
//             mirrorPlaneTransformation() const;

        // Edit

            //- Correct weighting factors for moving mesh.
            //  Updated for MeshObject.  HJ, 30/Aug/2010
            virtual bool movePoints() const;

            //- Update mesh topology using the morph engine
            //  Updated for MeshObject.  HJ, 30/Aug/2010
            virtual bool updateMesh(const mapPolyMesh&) const;

    // Interpolation functions

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        tmp<Field<Type> > interpolate
        (
            const polyPatch& patch,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        Type interpolate
        (
            const label pointID,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Get neighbour processor boundary face data
        template<class Type>
        void pointNgbProcBndFaceFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Get neighbour processor cell data
        template<class Type>
        void globalPointNgbProcBndFaceFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Get neighbour processor cell data
        template<class Type>
        void globalPointNgbProcCellFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Send and receive field data from neighbour processors
        template<class Type>
        tmp<FieldField<Field, Type> > procCellsFieldData
        (
            const Field<Type>& psi
        ) const;

        //- Send and receive field data from neighbour processors
        template<class Type>
        tmp<FieldField<Field, Type> > procBndFacesFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "leastSquaresVolPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
