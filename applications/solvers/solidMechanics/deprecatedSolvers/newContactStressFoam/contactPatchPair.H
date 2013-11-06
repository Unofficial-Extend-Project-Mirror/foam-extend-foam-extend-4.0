/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    contactPatchPair

Description
    A pair of surfaces in contact.
    
SourceFiles
    contactPatchPair.C

\*---------------------------------------------------------------------------*/

#ifndef contactPatchPair_H
#define contactPatchPair_H

#include "polyPatchID.H"
#include "dimensionedTypes.H"
#include "volFieldsFwd.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class contactProblem;

/*---------------------------------------------------------------------------*\
                       Class contactPatchPair Declaration
\*---------------------------------------------------------------------------*/

class contactPatchPair
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to contact problem
        const contactProblem& cp_;

        //- Master patch ID.  Gradient condition will be enforced
        const polyPatchID masterPatch_;

        //- Slave patch ID. Direction mixed condition will be enforced
        const polyPatchID slavePatch_;

        //- Friction coefficient
        dimensionedScalar frictionCoeff_;

        //- Contact tolerance
        const scalar contactTol_;

        //- Master patch interpolator
        primitivePatchInterpolation masterInterpolate_;

        //- Slave patch interpolator
        primitivePatchInterpolation slaveInterpolate_;

        //- Master to slave patch interpolator
        patchToPatchInterpolation masterToSlaveInterpolate_;

        //- Master to slave patch interpolator
        patchToPatchInterpolation slaveToMasterInterpolate_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        contactPatchPair(const contactPatchPair&);

        //- Disallow default bitwise assignment
        void operator=(const contactPatchPair&);



public:

    // Constructors

        //- Construct from components
        contactPatchPair
        (
            const word& name,
            const contactProblem& cp,
            const word& masterPatchName,
            const word& slavePatchName,
            const dimensionedScalar& frictionCoeff,
            const scalar contactTol,
            const intersection::algorithm alg = intersection::FULL_RAY,
            const intersection::direction dir = intersection::CONTACT_SPHERE
        );

        //- Construct from components
        contactPatchPair
        (
            const word& name,
            const contactProblem& cp,
            const dictionary& dict
        );


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return master patch ID
        const polyPatchID masterPatch() const
        {
            return masterPatch_;
        }

        //- Return slave patch ID
        const polyPatchID& slavePatch() const
        {
            return slavePatch_;
        }

        //- Return master touch fraction
        tmp<scalarField> masterTouchFraction() const;

        //- Return slave touch fraction
        tmp<scalarField> slaveTouchFraction() const;

        //- Correct contact data
        void correct
        (
            const FieldField<Field, vector>& curTraction,
            FieldField<Field, vector>& newTraction,
            FieldField<Field, vector>& refValue,
            FieldField<Field, scalar>& normalValueFraction
        );

        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
