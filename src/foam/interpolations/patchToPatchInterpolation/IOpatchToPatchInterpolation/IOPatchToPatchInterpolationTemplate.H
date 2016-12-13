/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
#include "PatchToPatchInterpolationTemplate.H"
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
#   include "IOPatchToPatchInterpolationTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
