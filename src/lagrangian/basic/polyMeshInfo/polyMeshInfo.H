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
    Foam::polyMeshInfo

Description


SourceFiles
    polyMeshInfo.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshInfo_H
#define polyMeshInfo_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class polyMeshInfo Declaration
\*---------------------------------------------------------------------------*/

class polyMeshInfo
{

    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Vector of valid directions in mesh
        //  defined according to the presence of empty patches
        mutable Vector<label> directions_;

        //- Number of number of valid geoemtric dimensions in the mesh
        label nGeometricD_;

        //- Number of number of valid solution dimensions in the mesh
        //  For wedge cases this includes the circumferential direction
        //  in case of swirl
        label nSolutionD_;

        //- Plane normal to the wedge patches, or between empty patches
        // for slab 2-D cases
        vector centreNormal_;

        //- Empty/unused component for 2-D cases
        label emptyComponent_;

        //- Geometric centre point
        vector centrePoint_;

        //- Number of empty patches
        label nEmpty_;

        //- Number of wedge patches
        label nWedge_;

        //- Axis of wedge
        vector wedgeAxis_;

        //- Angle of wedge
        scalar wedgeAngle_;


    // Private member functions

        //- Set the value of the mesh centre point
        void setCentrePoint();

        //- Determine wedge information
        void queryWedge();

        //- Determine mesh directions
        void queryDirections();


public:

    //- Runtime type information
    TypeName("polyMeshInfo");

    // Constructors

        //- Construct from components
        polyMeshInfo
        (
            const polyMesh& mesh
        );


    // Destructors

        virtual ~polyMeshInfo();


    // Member functions

        //- Access

            //- Return the number of valid geoemtric dimensions in the mesh
            label nGeometricD() const;

            //- Return the number of valid solution dimensions in the mesh
            label nSolutionD() const;

            //- Flag to indicate whether case is a 2-D wedge
            bool caseIs2dWedge() const;

            //- Flag to indicate whether case is a 2-D slab
            bool caseIs2dSlab() const;

            //- Flag to indicate whether case is 2-D
            bool caseIs2d() const;

            //- Flag to indicate whether case is 3-D
            bool caseIs3d() const;

            //- Return axis of wedge
            vector wedgeAxis() const;

            //- Return angle of wedge
            scalar wedgeAngle() const;

            //- Return the centre normal
            vector centreNormal() const;

            //- Return empty component
            label emptyComponent() const;

            //- Return centre point
            vector centrePoint() const;

            //- Set the empty component direction based on normal vector dir
            void setEmptyComponent(const vector& dir);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

