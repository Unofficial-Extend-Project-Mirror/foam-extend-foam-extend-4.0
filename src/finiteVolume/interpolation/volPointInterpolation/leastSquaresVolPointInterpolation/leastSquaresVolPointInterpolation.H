/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Interpolate volVectorField to mesh points using least squares best fit
    planes.

Please cite:
    Cardiff P, Karac A & Ivankovic A, A Large Strain Finite Volume Method for
    Orthotropic Bodies with General Material Orientations, Computer Methods
    in Applied Mechanics & Engineering, 2013,
    http://dx.doi.org/10.1016/j.cma.2013.09.008.

SourceFiles
    leastSquaresVolPointInterpolation.C

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVolPointInterpolation_H
#define leastSquaresVolPointInterpolation_H

#include "MeshObject.H"
#include "volFields.H"
#include "pointFields.H"
#include "simpleMatrix.H"

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

         //- mesh
         const fvMesh& mesh_;

        //- A matrix (in Ax=B) for each point
        //- A is the same for all components of phi
        //- A is a 4x4 for each point
        //List<scalarSquareMatrix> A_;

        //- B vector for each point
        //- each point has a 4 vectors
        //List<Field<vector> > B_;

    // Private member functions

        //- calculate A matrix for each point
        //- A is the same for all components of  phi
        void calcA(List<scalarSquareMatrix>& A) const;

        //- calc B source for each point
       void calcB
       (
           List<Field<vector> >& B,
           const GeometricField<vector, fvPatchField, volMesh>&
       ) const;

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

        //- Constructor given fvMesh.
        explicit leastSquaresVolPointInterpolation(const fvMesh&);


    // Destructor

        virtual ~leastSquaresVolPointInterpolation();


    // Member functions

        // Access


        // Edit

            virtual bool movePoints() const
            {
                return true;
            }

            virtual bool updateMesh(const mapPolyMesh&) const
            {
                return true;
            }


    // Interpolation functions

        //- Interpolate phi internal field and boundary field from volField
        //- to pointField using least squares fitted plane
        void interpolate
        (
            const GeometricField<vector, fvPatchField, volMesh>& vf,
            GeometricField<vector, pointPatchField, pointMesh>& pf
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
