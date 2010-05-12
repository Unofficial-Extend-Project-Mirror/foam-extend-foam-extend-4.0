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
    IOPatchToPatchInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

SourceFiles
    IOPatchToPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef IOPatchToPatchInterpolation_H
#define IOPatchToPatchInterpolation_H

#include "className.H"
#include "labelList.H"
#include "scalarField.H"
#include "pointField.H"
#include "FieldFields.H"
#include "faceList.H"
#include "intersection.H"
#include "PatchToPatchInterpolation.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class IOPatchToPatchInterpolationName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(IOPatchToPatchInterpolation);


/*---------------------------------------------------------------------------*\
                  Class IOPatchToPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class FromPatch, class ToPatch>
class IOPatchToPatchInterpolation
:
    public regIOobject,
    public PatchToPatchInterpolation<FromPatch, ToPatch>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        IOPatchToPatchInterpolation(const IOPatchToPatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const IOPatchToPatchInterpolation&);


public:

    // Constructors

        //- Construct from components
        IOPatchToPatchInterpolation
        (
            const IOobject& io,
            const FromPatch& fromPatch,
            const ToPatch& toPatch,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        );

    // Destructor

        ~IOPatchToPatchInterpolation();


    // Member Functions

        bool writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "IOPatchToPatchInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
