/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    ImmersedFace

Description
    Face intesection tools based on a distance function parameter
    Wet face has negative distance to the surface (ie it is below the surface)

SourceFiles
    ImmersedFace.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark
    Rewrite by Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef immersedFace_H
#define immersedFace_H

#include "polyMesh.H"
#include "immersedPoly.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ImmersedFace Declaration
\*---------------------------------------------------------------------------*/

template<class Distance>
class ImmersedFace
:
    public immersedPoly
{
    // Private data

        //- Reference to distance object
        const Distance& dist_;

        //- Point field with original face points and intersection points
        pointField facePointsAndIntersections_;

        //- Face that describes the submerged part of the original face
        face wetSubface_;

        //- Face that describes the dry part of the original face
        face drySubface_;

        //- All wet
        bool isAllWet_;

        //- All dry
        bool isAllDry_;


    // Private member functions

        //- Create points where the water surface intersects face edges
        void insertIntersectionPoints
        (
            const face& origFace,
            const scalarField& depth
        );

        //- Create wetSubface and drySubface
        void createSubfaces(const face& origFace, const scalarField& depth);


public:

    // Constructors

        //- Construct from components
        ImmersedFace
        (
            const label faceID,
            const polyMesh& mesh,
            const Distance& dist
        );


    //- Destructor - default


    // Public member functions

        // Access to primitive mesh data

            //- Return points
            inline const pointField& points() const
            {
                return facePointsAndIntersections_;
            }

            //- Return wet sub-face
            inline const face& wetSubFace() const
            {
                return wetSubface_;
            }

            //- Return dry sub-face
            inline const face& drySubFace() const
            {
                return drySubface_;
            }



        //- Is the face all wet?
        inline bool isAllWet() const
        {
            return isAllWet_;
        }

        //- Is the face all dry?
        inline bool isAllDry() const
        {
            return isAllDry_;
        }

        //- Return wet area magnitude
        inline scalar wetAreaMag() const
        {
            if (wetSubface_.empty())
            {
                // Wet face is empty if face is fully dry
                return 0;
            }
            else
            {
                return wetSubface_.mag(facePointsAndIntersections_);
            }
        }

        //- Return wet area centre
        inline point wetAreaCentre() const
        {
            if (wetSubface_.empty())
            {
                // Wet face is empty if face is fully dry
                return vector::zero;
            }
            else
            {
                return wetSubface_.centre(facePointsAndIntersections_);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "ImmersedFace.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
