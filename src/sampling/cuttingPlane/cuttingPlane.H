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
    Foam::cuttingPlane

Description
    Constructs plane through mesh.

    No attempt at resolving degenerate cases.

SourceFiles
    cuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPlane_H
#define cuttingPlane_H

#include "plane.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class cuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlane
:
    public plane
{
    // Private data

        //- List of cells cut by the plane
        labelList cutCells_;

        //- Intersection points
        pointField cutPoints_;

        //- Cut faces in terms of intersection points
        faceList cutFaces_;


    // Private Member Functions

        //- Determine cut cells, possibly restricted to a list of cells
        void calcCutCells
        (
            const primitiveMesh&,
            const scalarField& dotProducts,
            const labelList& cellIdLabels = labelList::null()
        );

        //- Determine intersection points (cutPoints_).
        labelList intersectEdges
        (
            const primitiveMesh&,
            const scalarField& dotProducts
        );

        //- Walk around circumference of cell starting from startEdgeI crossing
        //  only cut edges. Record cutPoint labels in faceVerts.
        static bool walkCell
        (
            const primitiveMesh&,
            const labelList& edgePoint,
            const label cellI,
            const label startEdgeI,
            DynamicList<label>& faceVerts
        );

        //- Determine cuts for all cut cells.
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const labelList& edgePoint
        );


protected:

    // Constructors

        //- Construct plane description without cutting
        cuttingPlane(const plane&);


    // Protected Member Functions

        //- recut mesh with existing planeDesc
        void reCut(const primitiveMesh&);

        //- recut mesh with existing planeDesc, restricted to a list of cells
        void reCut
        (
            const primitiveMesh&,
            const labelList& cellIdLabels
        );


public:

    // Constructors

        //- Construct from components: Mesh reference and plane
        cuttingPlane(const primitiveMesh&, const plane&);

        //- Construct from mesh reference and plane,
        //  restricted to a list of cells
        cuttingPlane
        (
            const primitiveMesh&,
            const plane&,
            const labelList& cellIdLabels
        );


    // Member Functions

        //- Return plane used
        const plane& planeDesc() const;

        //- Return pointField of cutting points
        const pointField& points() const;

        //- Return faceList of points in cells
        const faceList& faces() const;

        //- Return labelList of cut cells
        const labelList& cells() const;

        //- Return true or false to question: have any cells been cut?
        bool cut();

        //- Sample the cell field
        template<class Type>
        tmp<Field<Type> > sample(const Field<Type>&) const;

        template<class Type>
        tmp<Field<Type> > sample(const tmp<Field<Type> >&) const;


    // Member Operators

        void operator=(const cuttingPlane&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cuttingPlaneSample.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
