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
    Foam::surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolation_H
#define surfaceInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

        //- Discretisation schemes
        fvSchemes schemesDict_;

        //- Solver settings
        fvSolution solutionDict_;


        // Demand-driven data

            //- Central-differencing weighting factors
            mutable surfaceScalarField* weightingFactors_;

            //- Face-gradient difference factors
            mutable surfaceScalarField* differenceFactors_;

            //- Is mesh orthogonal
            mutable bool orthogonal_;

            //- Non-orthogonality correction vectors
            mutable surfaceVectorField* correctionVectors_;


    // Private member functions

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeCorrectionVectors() const;


protected:

    // Protected member functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and its debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    // Destructor

        ~surfaceInterpolation();


    // Member functions

        //- Return schemes
        const fvSchemes& schemesDict() const
        {
            return schemesDict_;
        }

        //- Return access to schemes
        fvSchemes& schemesDict()
        {
            return schemesDict_;
        }

        //- Return solver settings
        const fvSolution& solutionDict() const
        {
            return solutionDict_;
        }

        //- Return access to solver settings
        fvSolution& solutionDict()
        {
            return solutionDict_;
        }

        //- Return reference to weighting factors array
        const surfaceScalarField& weights() const;

        //- Return reference to difference factors array
        const surfaceScalarField& deltaCoeffs() const;

        //- Return whether mesh is orthogonal or not
        bool orthogonal() const;

        //- Return reference to non-orthogonality correction vectors array
        const surfaceVectorField& correctionVectors() const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
