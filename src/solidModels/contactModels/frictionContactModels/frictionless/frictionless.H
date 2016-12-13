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
    frictionless

Description
    frictionless contact friction less

SourceFiles
    frictionless.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef frictionless_H
#define frictionless_H

#include "frictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class frictionless Declaration
\*---------------------------------------------------------------------------*/

class frictionless
:
    public frictionContactModel
{
    // Private data

       //- slave friction disp - not used for frictionless
       vectorField slaveDisp_;

       //- slave friction traction - always zero for frictionless
       vectorField slaveTraction_;

       //- slave friction value fraction - always zero for frictionless
       symmTensorField slaveValueFrac_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        frictionless(const frictionless&);

        //- Disallow default bitwise assignment
        void operator=(const frictionless&);

public:

    //- Runtime type information
    TypeName("frictionless");

    // Static data members


    // Constructors

        //- Construct from dictionary
        frictionless
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID
        );


    // Destructor

        virtual ~frictionless();


    // Member Functions

        //- Correct the rheological model
        virtual void correct
        (
            const vectorField& slavePressure,
            const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
            const intersection::algorithm alg,
            const intersection::direction dir,
            const word interpolationMethod,
            const word fieldName,
            const Switch orthotropic,
            const word nonLinear,
            const vectorField& slaveFaceNormals
            )
        {}

        //- Return slave friction displacement
        virtual const vectorField& slaveDisp() const
        {
            return slaveDisp_;
        }

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slave friction value fraction
        virtual const symmTensorField& slaveValueFrac() const
        {
            return slaveValueFrac_;
        }

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
