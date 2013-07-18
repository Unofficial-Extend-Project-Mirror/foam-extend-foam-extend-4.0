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

Author
    philip.cardiff@ucd.ie
    philip.cardiff@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef contactPatchPair_H
#define contactPatchPair_H

#include "polyPatchID.H"
#include "dimensionedTypes.H"
#include "volFieldsFwd.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "ggiInterpolation.H"
#include "FieldFields.H"
#include "pointFields.H"

#include "fvMesh.H"

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

        //- contact on or off
        Switch contactActive_;

        //- Name
        const word name_;

        //- Reference to contact problem
        contactProblem& cp_;

        //- patch to patch algorithm
        const intersection::algorithm alg_;

        //- patch to patch direction
        const intersection::direction dir_;

        //- Master patch ID.
        const polyPatchID masterPatch_;

        //- Slave patch ID.
        const polyPatchID slavePatch_;

        //- Contact gap tolerance
        const scalar contactGapTol_;

        //- penalty scale factor
        scalar penaltyScaleFactor_;

        //- penalty factor
        scalar returnScaleFactor_;

        //contact iteration number
        int contactIterNum_;

    // file that contains contact information
        OFstream* contactFilePtr_;

        // slave face zone name - needed to define zoneToZoneInterp
        word slaveFaceZoneName_;

        // master face zone name - needed to define zoneToZoneInterp
        word masterFaceZoneName_;

        // slave point zone name - needed to define the size of global point fields
        word slavePointZoneName_;

        // accumulated slave point force
        vectorField totalSlavePointForce_;
  //IOField<vector> totalSlavePointForce_;

        // current slave point penetration
        scalarField slavePointPenetration_;

        // current master point penetration
        scalarField masterPointPenetration_;

        // current global slave point penetration
        scalarField globalSlavePointPenetration_;

        // current global master point penetration
        scalarField globalMasterPointPenetration_;

        // accumulated slave point force
        vectorField oldTotalSlavePointForce_;

        // old minimum slave point penetration
        scalar oldMinSlavePointPenetration_;

        // penetration limit - for when the master surface surrounds the slave
        scalar penetrationLimit_;

        //- whether the master surface is considered rigid
        Switch rigidMaster_;

        //- interpolation method for passing traction from slave to master
        word interpolationMethod_;

       //- zoneToZoneInterpolation for the calculation of pointDistanceToIntersection
       zoneToZoneInterpolation faceZoneMasterToSlaveInterpolator_;

       //- zoneToZone or ggiZone for interpolation of traction from slave to master
       zoneToZoneInterpolation* faceZoneSlaveToMasterInterpolatorPtr_;
       ggiZoneInterpolation* faceZoneGgiInterpolatorPtr_;

       //- interpolation from slave vertices to slave face centres
       primitivePatchInterpolation slaveInterpolator_;


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
            contactProblem& cp,
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

        //- Return master point gap
        tmp<scalarField> masterGapPoints() const;

        //- Return slave point gap
        tmp<scalarField> slaveGapPoints() const;

        //- Return master total point force
        tmp<vectorField> masterPointForce() const;

        //- Return slave total point force
        tmp<vectorField> slavePointForce() const;

        //- Return master contact pressure for visualisation
        tmp<scalarField> masterContactPressure() const;

        //- Return slave contact pressure for visualisation
        tmp<scalarField> slaveContactPressure() const;

        //- Correct contact data
        void correct();

        //- if contact correction on or off
        Switch contactActive()
        {
             return contactActive_;
        };

        //- set contactActive_ to on or off
        void setContactActive(bool trueOrFalse)
        {
            contactActive_ = trueOrFalse;
        };

        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
