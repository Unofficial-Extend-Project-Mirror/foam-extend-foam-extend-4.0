/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
    labelHashSet* setPtr = nullptr
);

//- Check for negative cell volumes
bool checkCellVolumes
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = nullptr
);

//- Check for negative face areas
bool checkFaceAreas
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minFaceArea = VSMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);


//- Check quality of tetrahedra
bool checkTetQuality
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minTetQuality = VSMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);


//- Check minimum face twist
bool checkMinTwist
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minTwist = VSMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);


//- Check the area of internal faces versus boundary faces
bool checkCellDeterminant
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnDet = 1e-3,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check volume ratio
void checkMinVolRatio
(
    const polyMeshGen&,
    scalarField&,
    const boolList* changedFacePtr = nullptr
);

bool checkMinVolRatio
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnVolRatio = 0.01,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face triangle twist
bool checkTriangleTwist
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minTwist = VSMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check for negative part tetrahedra
//- Cells are decomposed into tetrahedra formed by
//- the cell centre, face centre and the edge vertices
bool checkCellPartTetrahedra
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minPartTet = VSMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check for non-orthogonality
void checkFaceDotProduct
(
    const polyMeshGen&,
    scalarField&,
    const boolList* changedFacePtr = nullptr
);

bool checkFaceDotProduct
(
    const polyMeshGen&,
    const bool report = false,
    const scalar nonOrthWarn = 70.0,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face pyramid volume
bool checkFacePyramids
(
    const polyMeshGen&,
    const bool report = false,
    const scalar minPyrVol = -SMALL,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face skewness
void checkFaceSkewness
(
    const polyMeshGen&,
    scalarField&,
    const boolList* changedFacePtr = nullptr
);

bool checkFaceSkewness
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnSkew = 4.0,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face uniformity
void checkFaceUniformity
(
    const polyMeshGen&,
    scalarField&,
    const boolList* changedFacePtr = nullptr
);

bool checkFaceUniformity
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnUniform = 0.1,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- check volume difference of neighbouring cells
void checkVolumeUniformity
(
    const polyMeshGen&,
    scalarField&,
    const boolList* changedFacePtr = nullptr
);

bool checkVolumeUniformity
(
    const polyMeshGen&,
    const bool report = false,
    const scalar warnUniform = 0.1,
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face angles
bool checkFaceAngles
(
    const polyMeshGen&,
    const bool report = false,
    const scalar maxDeg = 10,    // In degrees
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

//- Check face warpage: decompose face and check ratio between
//  magnitude of sum of triangle areas and sum of magnitude of
//  triangle areas.
bool checkFaceFlatness
(
    const polyMeshGen&,
    const bool report,
    const scalar warnFlatness,  // When to include in set.
    labelHashSet* setPtr = nullptr,
    const boolList* changedFacePtr = nullptr
);

// Checks using topology only

//- Check for unused points
bool checkPoints
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = nullptr
);

//- Check face ordering
bool checkUpperTriangular
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = nullptr
);

//- Check cell zip-up
bool checkCellsZipUp
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = nullptr
);

//- Check uniqueness of face vertices
bool checkFaceVertices
(
    const polyMeshGen&,
    const bool report = false,
    labelHashSet* setPtr = nullptr
);

//- Check mesh topology for correctness. Returns false for no error.
bool checkTopology(const polyMeshGen&, const bool report = false);

//- Check mesh geometry (& implicitly topology) for correctness.
//  Returns false for no error.
bool checkGeometry(const polyMeshGen&, const bool report = false);

//- Check mesh for correctness. Returns false for no error.
bool checkMesh(const polyMeshGen&, const bool report = false);

//- Read the user defined mesh quality settings
label findBadFacesAdditionalChecks
(
    const polyMeshGen& mesh,
    const bool report,
    labelHashSet& badFaces,
    const boolList* activeFacePtr = nullptr
);

//- checks for bad faces making the mesh unusable
//- checks for negative pyramids and zero area faces
label findBadFacesRelaxed
(
    const polyMeshGen&,
    labelHashSet& badFaces,
    const bool report = false,
    const boolList* activeFacePtr = nullptr
);

//- check mesh for bad faces which make the mesh invalid
//- checks for negative pyramids, face flatness,
//- zero area faces and negative tetrahedra
label findBadFaces
(
    const polyMeshGen&,
    labelHashSet& badFaces,
    const bool report = false,
    const boolList* activeFacePtr = nullptr
);

//- check mesh for faces that may cause problems to the solver
//- checks for non-orthogonality and skewness
label findLowQualityFaces
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool report = false,
    const boolList* activeFacePtr = nullptr
);

//- checks the mesh and selects the faces with worst quality
label findWorstQualityFaces
(
    const polyMeshGen& mesh,
    labelHashSet& badFaces,
    const bool report = false,
    const boolList* activeFacePtr = nullptr,
    const scalar relativeThreshold = 0.1
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
