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
    Foam::meshEdgeCuts

Description
    Container for cuts of edges of mesh. Cuts either exactly through existing 
    vertices or through edges.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshEdgeCuts_H
#define meshEdgeCuts_H

#include "labelList.H"
#include "labelHashSet.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class meshEdgeCuts Declaration
\*---------------------------------------------------------------------------*/

class meshEdgeCuts
{
    // Private data

        const primitiveMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshEdgeCuts(const meshEdgeCuts&);

        //- Disallow default bitwise assignment
        void operator=(const meshEdgeCuts&);


protected:

        //- List of cells containing the cuts
        labelList cells_;

        //- Points exactly cut by cuts
        labelList meshVerts_;

        //- List of edge labels cut
        labelList meshEdges_;

        //- Positions on edges
        scalarField meshEdgeWeights_;


    // Protected Static Functions

        //- Mark element in a hashSet
        static void mark(const label elem, labelHashSet& markedElems);

        //- Mark list of elements in a hashSet
        static void mark(const labelList& elems, labelHashSet& markedElems);

        //- Return true and set weight if linear interpolation between
        //  val0 and val1 crosses isoVal. weight=1 if isoVal==val1
        static bool crosses
        (
            const scalar isoVal,
            const scalar val0,
            const scalar val1,
            scalar& weight
        );

public:

    // Constructors

        //- Construct from components
        meshEdgeCuts
        (
            const primitiveMesh& mesh,
            const labelList& cells,
            const labelList& meshVerts,
            const labelList& meshEdges,
            const scalarField& meshEdgeWeights
        );


    // Member Functions

        const primitiveMesh& mesh() const
        {
            return mesh_;
        }

        const labelList& cells() const
        {
            return cells_;
        }

        const labelList& meshVerts() const
        {
            return meshVerts_;
        }

        const labelList& meshEdges() const
        {
            return meshEdges_;
        }

        const scalarField& meshEdgeWeights() const
        {
            return meshEdgeWeights_;
        }

        label size() const
        {
            return meshVerts_.size() + meshEdges_.size();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
