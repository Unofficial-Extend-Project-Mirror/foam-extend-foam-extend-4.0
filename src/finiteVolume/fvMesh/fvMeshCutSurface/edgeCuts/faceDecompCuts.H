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
    Foam::faceDecompCuts

Description
    Container for cuts of edges of (implicit) tet decomposition. Used to
    collect data for meshCut. As much as possible cuts are defined using
    mesh information:

    - cut (exactly) through mesh vertex
    -   ,,              ,,       cell centre
    -   ,,              ,,       face centre

    - cut through mesh edge. Both edge label and position on edge given.

    - cut through tet pyramidEdge (edge between vertex and cell centre).
      Edge and position on edge given.
    - cut through tet faceEdge (edge between vertex and face centre).
      Edge and position on edge given.
    - cut through tet centreEdge (edge between face centre and cell centre).
      Edge and position on edge given.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef faceDecompCuts_H
#define faceDecompCuts_H

#include "meshEdgeCuts.H"
#include "pyramidEdge.H"
#include "faceEdge.H"
#include "centreEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class faceDecompCuts Declaration
\*---------------------------------------------------------------------------*/

class faceDecompCuts
:
    public meshEdgeCuts
{

protected:

        //- List of faces whose centre is exactly cut
        labelList meshFaceCentres_;

        //- List of cells whose centre is exactly cut
        labelList meshCellCentres_;

        //- List of pyramid edge descriptions
        List<pyramidEdge> pyrEdges_;

        //- For each pyramid edge description the position of the cut
        scalarField pyrEdgeWeights_;

        //- List of face edge descriptions
        List<faceEdge> faceEdges_;

        //- For each face edge description the position of the cut
        scalarField faceEdgeWeights_;

        //- List of centre edge descriptions
        List<centreEdge> centreEdges_;

        //- For each centre edge description the position of the cut
        scalarField centreEdgeWeights_;

public:

    // Constructors

        //- Construct from components
        faceDecompCuts
        (
            const primitiveMesh& mesh,
            const labelList& cells,

            const labelList& meshVerts,
            const labelList& meshFaceCentres,
            const labelList& meshCellCentres,

            const labelList& meshEdges,
            const scalarField& meshEdgeWeights,

            const List<pyramidEdge>& pyrEdges,
            const scalarField& pyrEdgeWeights,

            const List<faceEdge>& faceEdges,
            const scalarField& faceEdgeWeights,

            const List<centreEdge>& centreEdges,
            const scalarField& centreEdgeWeights
        )
        :
            meshEdgeCuts(mesh, cells, meshVerts, meshEdges, meshEdgeWeights),
            meshFaceCentres_(meshFaceCentres),
            meshCellCentres_(meshCellCentres),
            pyrEdges_(pyrEdges),
            pyrEdgeWeights_(pyrEdgeWeights),
            faceEdges_(faceEdges),
            faceEdgeWeights_(faceEdgeWeights),
            centreEdges_(centreEdges),
            centreEdgeWeights_(centreEdgeWeights)
        {}


    // Member Functions

        const labelList& meshFaceCentres() const
        {
            return meshFaceCentres_;
        }

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

        const List<centreEdge>& centreEdges() const
        {
            return centreEdges_;
        }

        const scalarField& centreEdgeWeights() const
        {
            return centreEdgeWeights_;
        }

        const List<faceEdge>& faceEdges() const
        {
            return faceEdges_;
        }

        const scalarField& faceEdgeWeights() const
        {
            return faceEdgeWeights_;
        }

        label size() const
        {
            return
                meshEdgeCuts::size()
              + meshFaceCentres_.size() + meshCellCentres_.size()
              + pyrEdges_.size() + centreEdges_.size() + faceEdges_.size();       
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
