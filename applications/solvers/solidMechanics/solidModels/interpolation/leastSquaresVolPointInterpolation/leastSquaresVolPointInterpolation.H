/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::leastSquaresVolPointInterpolation

Description
    Interpolate volVectorField to mesh points using least squares best fit
    planes

SourceFiles
    leastSquaresVolPointInterpolation.C

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVolPointInterpolation_H
#define leastSquaresVolPointInterpolation_H

#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class leastSquaresVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/
class leastSquaresVolPointInterpolation
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
       void calcB(List<Field<vector> >& B, const GeometricField<vector, fvPatchField, volMesh>&) const;

        //- Disallow default bitwise copy construct
        leastSquaresVolPointInterpolation(const leastSquaresVolPointInterpolation&);

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
