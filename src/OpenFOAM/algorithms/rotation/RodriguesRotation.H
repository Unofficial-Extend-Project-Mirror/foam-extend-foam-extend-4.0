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
    RodriguesRotation

Description
    Compute rotation tensor base on Rodrigues' Rotation Formula
    http://mathworld.wolfram.com/RodriguesRotationFormula.html
    http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

Author
    Martin Beaudoin, Hydro-Quebec, (2008, 2011)

SourceFiles
    RodriguesRotation.C

\*---------------------------------------------------------------------------*/

#ifndef RodriguesRotation_H
#define RodriguesRotation_H

#include "tmp.H"
#include "vector.H"
#include "tensor.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Compute rotation tensor using Rodrigues' Rotation Formula
// Formulation taken from:
// http://mathworld.wolfram.com/RodriguesRotationFormula.html
//
// Requires: rotation axis + rotation angle
tensor RodriguesRotation
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees = true
);

// Compute rotation tensor using Rodrigues' Rotation Formula
// Slightly different formulation taken from:
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + rotation angle
tensor RodriguesRotation2
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees = true
);

// Compute rotation tensor using Rodrigues' Rotation Formula
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + initial vector + final vector

tensor RodriguesRotation
(
    const vector& rotationAxis,
    const vector& vi,
    const vector& vf
);

// Compute rotation tensorField using Rodrigues' Rotation Formula
// http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
// Requires: rotation axis + initial vectorField + final vectorField
Foam::tensorField RodriguesRotation
(
    const vector& rotationAxis,
    const vectorField& vi,
    const vectorField& vf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
