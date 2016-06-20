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
    crackerMovingFvMesh

Description

SourceFiles
    crackerMovingFvMesh.C

Author
    Zeljko Tukovic, FSB Zagreb  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef crackerMovingFvMesh_H
#define crackerMovingFvMesh_H

#include "crackerFvMesh.H"
#include "dimensionedTypes.H"
#include "laplaceTetMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class laplaceTetMotionSolver;

/*---------------------------------------------------------------------------*\
                   Class crackFvMesh Declaration
\*---------------------------------------------------------------------------*/

class crackerMovingFvMesh
:
    public crackerFvMesh
{
    // Private data

//         //- Moving mesh dictionary
//         dictionary dynamicMeshCoeffs_;

        //- Moving mesh solver
        autoPtr<laplaceTetMotionSolver> motionPtr_;

        //- Notch patch name
        word notchPatchName_;

        //- Notch patch index
        label notchPatchIndex_;

        //-
        label crackTipEdgeIndex_;

//         //-
//         vector crackTipEdgeNormal_;

        //-
        label maxCosFaceIndex_;



        //- Spring analogy mesh motion
        Switch springMeshMotion_;

        //-
        scalar deltaT_;

        //-
        scalar dpTol_;

        //-
        label maxIters_;

        //- List of constrained points
        DynamicList<label> fixedPoints_;

        //- Old points
        vectorField oldPoints_;

    // Private Member Functions

        //- Mesh size function
        scalar sizeFunction(const vector& p) const;

        //- Disallow default bitwise copy construct
        crackerMovingFvMesh(const crackerMovingFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const crackerMovingFvMesh&);

public:

    //- Runtime type information
    TypeName("crackerMovingFvMesh");


    // Constructors

        //- Construct from IOobject
        crackerMovingFvMesh(const IOobject& io);


    // Destructor

        ~crackerMovingFvMesh();


    // Member Functions


        //- Return notch patch index
        label notchPatchIndex() const
        {
            return notchPatchIndex_;
        }

        //- Update crack tip
        bool updateCrackTip();

        //- Smooth mesh (internal points)
        void smooth(vectorField& points);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update the mesh for both mesh motion and topology change
        virtual bool smoothMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
