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

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

\*---------------------------------------------------------------------------*/

#include "IOPatchToPatchInterpolation.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class FromPatch, class ToPatch>
IOPatchToPatchInterpolation<FromPatch, ToPatch>::IOPatchToPatchInterpolation
(
    const IOobject& io,
    const FromPatch& fromPatch,
    const ToPatch& toPatch,
    intersection::algorithm alg,
    const intersection::direction dir
)
:
    regIOobject(io),
    PatchToPatchInterpolation<FromPatch, ToPatch>(fromPatch, toPatch, alg, dir)
{
    if (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    {
        Istream& is = readStream(typeName);

        labelList* paPtr = new labelList(is);
        FieldField<Field, scalar>* pwPtr = new FieldField<Field, scalar>(is);
        scalarField* pdPtr = new scalarField(is);

        labelList* faPtr = new labelList(is);
        FieldField<Field, scalar>* fwPtr = new FieldField<Field, scalar>(is);
        scalarField* fdPtr = new scalarField(is);
        Info << "Setting weights from file" << endl;
        this->setWeights(paPtr, pwPtr, pdPtr, faPtr, fwPtr, fdPtr);
    }
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
IOPatchToPatchInterpolation<FromPatch, ToPatch>::~IOPatchToPatchInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
bool IOPatchToPatchInterpolation<FromPatch, ToPatch>::writeData
(
    Ostream& os
) const
{
    os << this->pointAddr() << nl;
    os << this->pointWeights() << nl;
    os << this->pointDistanceToIntersection() << nl;
    os << this->faceAddr() << nl;
    os << this->faceWeights() << nl;
    os << this->faceDistanceToIntersection() << nl;
    return os.good();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
