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
    ImmersedCell

Description
    ImmersedCell is created from a mesh cell with a distance function.
    Signed distance is provided to indicate the cell cut into the
    "wet part" (with negative distance) and "dry part" (with positive distance).
    Based on the data, the cell is cut by the distance = 0 plane/face into
    two sub-cells and data for the wet cell part is provided.

Implementation
    The class builds a primitiveMesh with the local point, face and cell list
    which is then manipulated as needed.

SourceFiles
    ImmersedCell.C

Author
    Dominik Christ, Wikki Ltd.
    Based on waves2Foam by
    Niels Gjoel Jacobsen, Technical University of Denmark
    Rewrite and debugging Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef ImmersedCell_H
#define ImmersedCell_H

#include "immersedPoly.H"
#include "polyMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ImmersedCell Declaration
\*---------------------------------------------------------------------------*/

template<class Distance>
class ImmersedCell
:
    public immersedPoly,
    public primitiveMesh
{
    // Private data

        //- Global ID of the cell
        label cellID_;

        //- Local reference to mesh
        const polyMesh& mesh_;

        //- Reference to distance object
        const Distance& dist_;

        //- Point field with original cell points and intersection points
        pointField cellPointsAndIntersections_;

        //- Minimum edge length in cell
        scalar absTol_;

        //- All wet
        bool isAllWet_;

        //- All dry
        bool isAllDry_;


        // Primitive mesh data for two cell mesh

            // Mesh points
            pointField points_;

            // Mesh faces
            faceList faces_;

            // Face owners
            labelList faceOwner_;

            // Face neigbour
            labelList faceNeighbour_;

        //- Depth to immersed surface for cell points and intersections
        //  Note: positive depth is dry
        scalarField depth_;
    

    // Private member functions

        //- Copied from class geomCellLooper
        void getBase(const vector& n, vector& e0, vector& e1) const;

        //- Create points where the water surface intersects face edges
        //  Return number of inserted points
        //  and if so, insert additional points into points_ and adjust depth
        label insertIntersectionPoints();

        //- Make an internal face out of the intersection points
        face createInternalFace(const label nIntersections) const;


public:

    // Constructors

        //- Construct from components
        ImmersedCell
        (
            const label cellID,
            const polyMesh& mesh,
            const Distance& dist
        );


    //- Destructor - default


    // Public member functions

        // Access to primitive mesh data

            //- Return mesh points
            virtual const pointField& points() const
            {
                return points_;
            }

            //- Return faces
            virtual const faceList& faces() const
            {
                return faces_;
            }


            //- Face face-owner addresing
            virtual const labelList& faceOwner() const
            {
                return faceOwner_;
            }

            //- Face face-neighbour addressing
            virtual const labelList& faceNeighbour() const
            {
                return faceNeighbour_;
            }

            //- Return old points for mesh motion
            virtual const pointField& oldPoints() const
            {
                return points_;
            }


        //- Is the cell all wet?
        inline bool isAllWet() const
        {
            return isAllWet_;
        }

        //- Is the cell all dry?
        inline bool isAllDry() const
        {
            return isAllDry_;
        }

        //- Return center of the cell's wet part
        inline point wetVolumeCentre() const
        {
            if (isAllWet_ || isAllDry_)
            {
                return mesh_.cellCentres()[cellID_];
            }
            else
            {
                return this->cellCentres()[WET];
            }
        }

        //- Return volume of the cell's wet part
        inline scalar wetVolume() const
        {
            if (isAllWet_) // wet face is empty if face is fully dry
            {
                return mesh_.cellVolumes()[cellID_];
            }
            else if (isAllDry_)
            {
                return 0;
            }
            else
            {
                return this->cellVolumes()[WET];
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "ImmersedCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
