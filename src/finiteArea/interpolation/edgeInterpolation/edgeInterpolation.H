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
    edgeInterpolation

Description
    Cell to edge interpolation scheme. Included in faMesh.

SourceFiles
    edgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef edgeInterpolation_H
#define edgeInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "faSchemes.H"
#include "faSolution.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class edgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

class edgeInterpolation
:
    public faSchemes,
    public faSolution
{
    // Private data

        // Reference to faMesh
        const faMesh& faMesh_;

        //- Geodesic distance between centroides of neighbour finite areas
        mutable edgeScalarField* lPN_;

        //- Cenral-differencing weighting factors
        mutable edgeScalarField* weightingFactors_;

        //- Face-gradient difference factors
        mutable edgeScalarField* differenceFactors_;

        //- Is mesh orthogonal
        mutable bool orthogonal_;

        //- Non-orthogonality correction vectors
        mutable edgeVectorField* correctionVectors_;

        //- Is mesh skew
        mutable bool skew_;

        //- Skew correction vectors
        mutable edgeVectorField* skewCorrectionVectors_;


    // Private member functions

        //- Construct geodesic distance between P and N
        void makeLPN() const;

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeCorrectionVectors() const;

        //- Construct skewness correction vectors
        void makeSkewCorrectionVectors() const;

//         //- Construct Least-squares gradient vectors
//         void makeLeastSquareVectors() const;


protected:

    // Protected member functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and it's debug switch
    ClassName("edgeInterpolation");


    // Constructors

        //- Construct given an faMesh
        edgeInterpolation(const faMesh&);


    // Destructor

        ~edgeInterpolation();


    // Member functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return faMesh_;
        }

        //- Return reference to PN geodesic distance
        const edgeScalarField& lPN() const;

        //- Return reference to weighting factors array
        const edgeScalarField& weights() const;

        //- Return reference to difference factors array
        const edgeScalarField& deltaCoeffs() const;

        //- Return whether mesh is orthogonal or not
        bool orthogonal() const;

        //- Return reference to non-orthogonality correction vectors array
        const edgeVectorField& correctionVectors() const;

        //- Return whether mesh is skew or not
        bool skew() const;

        //- Return reference to skew vectors array
        const edgeVectorField& skewCorrectionVectors() const;


        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
