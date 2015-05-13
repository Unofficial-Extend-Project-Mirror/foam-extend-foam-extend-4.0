/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    polyMeshGenChecks

Description
    A set of functions used for mesh checking mesh quality

SourceFiles
    polyMeshGenAddressingChecks.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenChecks_H
#define polyMeshGenChecks_H

#include "polyMeshGen.H"
#include "boolList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Namespace polyMeshGenChecks functions Declaration
\*---------------------------------------------------------------------------*/

namespace polyMeshGenChecks
{
    
// Check geometry
//- Check boundary closedness
bool checkClosedBoundary(const polyMeshGen&, const bool report = false);

//- Check cells for closedness
bool checkClosedCells
(
    const polyMeshGen&,
    const bool report = false,
    const scalar aspectWarn = 1000,
    labelHashSet* setPtr = NULL
);

//- Check for negative cell volumes
bool checkCellVolumes
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = NULL
);

//- Check for negative face areas
bool checkFaceAreas
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minFaceArea = VSMALL,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check for negative part tetrahedra
//- Cells are decomposed into tetrahedra formed by
//- the cell centre, face centre and the edge vertices
bool checkCellPartTetrahedra
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minPartTet = VSMALL,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check for non-orthogonality
bool checkFaceDotProduct
(
    const polyMeshGen&,
    const bool report = false,
    const scalar nonOrthWarn = 70.0,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check face pyramid volume
bool checkFacePyramids
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minPyrVol = -SMALL,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check face skewness
bool checkFaceSkewness
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnSkew = 4.0,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check face uniformity
bool checkFaceUniformity
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnUniform = 0.1,
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check face angles
bool checkFaceAngles
(
    const polyMeshGen&,
    const bool report = false,
    const scalar maxDeg = 10,    // In degrees
    labelHashSet* setPtr = NULL,
    const boolList* changedFacePtr = NULL
);

//- Check face warpage: decompose face and check ratio between
//  magnitude of sum of triangle areas and sum of magnitude of
//  triangle areas.
bool checkFaceFlatness
(
    const polyMeshGen&,
    const bool report,
    const scalar warnFlatness,  // When to include in set.
    labelHashSet* setPtr,
    const boolList* changedFacePtr = NULL
);

// Checks using topology only

//- Check for unused points
bool checkPoints
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = NULL
);

//- Check face ordering
bool checkUpperTriangular
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = NULL
);

//- Check cell zip-up
bool checkCellsZipUp
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = NULL
);

//- Check uniqueness of face vertices
bool checkFaceVertices
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = NULL
);

//- Check mesh topology for correctness. Returns false for no error.
bool checkTopology(const polyMeshGen&, const bool report = false);

//- Check mesh geometry (& implicitly topology) for correctness. 
//  Returns false for no error.
bool checkGeometry(const polyMeshGen&, const bool report = false);

//- Check mesh for correctness. Returns false for no error.
bool checkMesh(const polyMeshGen&, const bool report = false);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
