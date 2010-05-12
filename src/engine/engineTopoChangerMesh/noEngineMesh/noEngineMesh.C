// The FOAM Project // File: noEngineMesh.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   noEngineMesh
   \\  /           | Family: engine
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR

-------------------------------------------------------------------------------
*/

#include "noEngineMesh.H"
#include "engineTime.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "volPointInterpolation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noEngineMesh, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, noEngineMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::noEngineMesh::noEngineMesh
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noEngineMesh::~noEngineMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
