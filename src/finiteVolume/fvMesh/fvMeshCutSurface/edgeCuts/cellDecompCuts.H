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
    Foam::cellDecompCuts

Description
    Container for cuts of edges of (implicit) tet decomposition.
    Used to collect data for meshCut.

    As much as possible, cuts are defined using mesh information:
      - cut (exactly) through mesh vertex
      - cut (exactly) through cell centre

      - cut through mesh edge. Both edge label and position on edge given.

      - cut through tet pyramidEdge (edge between vertex and cell centre).
        Edge and position on edge given.

      - cut through diagonalEdge (edge between vertices of a face)
        Edge and position on edge given.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef cellDecompCuts_H
#define cellDecompCuts_H

#include "meshEdgeCuts.H"
#include "pyramidEdge.H"
#include "diagonalEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class cellDecompCuts Declaration
\*---------------------------------------------------------------------------*/

class cellDecompCuts
:
    public meshEdgeCuts
{

protected:

        labelList meshCellCentres_;

        List<pyramidEdge> pyrEdges_;
        scalarField pyrEdgeWeights_;

        List<diagonalEdge> diagEdges_;
        scalarField diagEdgeWeights_;

    // Private Member Functions


public:

    // Constructors

        //- Construct from components
        cellDecompCuts
        (
            const primitiveMesh& mesh,
            const labelList& cells,

            const labelList& meshVerts,
            const labelList& meshCellCentres,

            const labelList& meshEdges,
            const scalarField& meshEdgeWeights,

            const List<pyramidEdge>& pyrEdges,
            const scalarField& pyrEdgeWeights,

            const List<diagonalEdge>& diagEdges,
            const scalarField& diagEdgeWeights
        )
        :
            meshEdgeCuts(mesh, cells, meshVerts, meshEdges, meshEdgeWeights),
            meshCellCentres_(meshCellCentres),
            pyrEdges_(pyrEdges),
            pyrEdgeWeights_(pyrEdgeWeights),
            diagEdges_(diagEdges),
            diagEdgeWeights_(diagEdgeWeights)
        {}


    // Member Functions

        const labelList& meshCellCentres() const
        {
            return meshCellCentres_;
        }

        const List<pyramidEdge>& pyrEdges() const
        {
            return pyrEdges_;
        }

        const scalarField& pyrEdgeWeights() const
        {
            return pyrEdgeWeights_;
        }

        const List<diagonalEdge>& diagEdges() const
        {
            return diagEdges_;
        }

        const scalarField& diagEdgeWeights() const
        {
            return diagEdgeWeights_;
        }

        label size() const
        {
            return
                meshEdgeCuts::size() + meshCellCentres_.size()
              + pyrEdges_.size() + diagEdges_.size();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
