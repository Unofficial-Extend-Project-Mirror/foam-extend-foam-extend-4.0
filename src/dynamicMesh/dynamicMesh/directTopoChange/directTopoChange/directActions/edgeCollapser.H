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
    Foam::edgeCollapser

Description
    Does directTopoChanges to remove edges. Can remove faces due to edge
    collapse but can not remove cells due to face removal!
    Also removes unused points.

SourceFiles
    edgeCollapser.C

\*---------------------------------------------------------------------------*/

#ifndef edgeCollapser_H
#define edgeCollapser_H

#include "labelList.H"
#include "dynamicLabelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class directTopoChange;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class edgeCollapser Declaration
\*---------------------------------------------------------------------------*/

class edgeCollapser
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- For every point -1 or region number
        labelList pointRegion_;

        //- -1 or master vertex for region number
        dynamicLabelList pointRegionMaster_;

        //- Stack of free region numbers. Corresponds to -1 in pointRegionMaster
        SLList<label> freeRegions_;


    // Static Functions

        //- Find val in list. Search starts at start, continues to size-1.
        static label findIndex
        (
            const labelList&,
            const label start,
            const label size,
            const label val
        );


    // Private Member Functions

        //- Determine points connected through edgesToRemove_.
        //  Note: Only routine that uses edgesToRemove!
        label changePointRegion
        (
            const label pointI,
            const label oldRegion,
            const label newRegion
        );

        //- Whether point is master of region or has been removed
        bool pointRemoved(const label) const;

        //- Renumber f with new vertices. Removes duplicates.
        void filterFace(const label faceI, face&) const;

        //- Some debugging printing
        void printRegions() const;

        //- Collapse list of edges. Tries to find master to collapse to.
        void collapseEdges(const labelList& edgeLabels);

        //- Disallow default bitwise copy construct
        edgeCollapser(const edgeCollapser&);

        //- Disallow default bitwise assignment
        void operator=(const edgeCollapser&);


public:

    //- Runtime type information
    ClassName("edgeCollapser");


    // Constructors

        //- Construct from mesh.
        edgeCollapser(const polyMesh& mesh);


    // Member Functions

        // Access

            //- For every point the region it belongs to or -1.
            const labelList& pointRegion() const
            {
                return pointRegion_;
            }

            //- For every region the master (i.e. the point the region will
            //  be replaced by)
            const dynamicLabelList& pointRegionMaster() const
            {
                return pointRegionMaster_;
            }

            //- Check that edge is not marked for anything
            bool unaffectedEdge(const label edgeI) const;


        // Edit

            //- Set edge to collapse and point to collapse it to.
            //  Return true if collapse is valid.
            //  (always true at the moment)
            bool collapseEdge(const label edgeI, const label master);

            //- Play commands into directTopoChange to create mesh. Return true
            //  if anything changed.
            bool setRefinement(directTopoChange&);

            void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
