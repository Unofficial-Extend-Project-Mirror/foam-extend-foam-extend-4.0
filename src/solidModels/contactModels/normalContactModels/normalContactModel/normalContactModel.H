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
    normalContactModel

Description
    model for the behaviour of the normal contact

SourceFiles
    normalContactModel.C
    newNormalContactModel.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef normalContactModel_H
#define normalContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"
#include "patchToPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "ggiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class normalContactModel Declaration
\*---------------------------------------------------------------------------*/

class normalContactModel
{
    // Private data

        //- Name
        const word name_;

        //- Reference to primitivePatch
        const fvPatch& patch_;

        //- Master and slave patch IDs
        const label masterPatchID_;
        const label slavePatchID_;

        //- Master and slave face zone ID
        const label masterFaceZoneID_;
        const label slaveFaceZoneID_;

        //- contact point gap
        scalarField slaveContactPointGap_;

       //- zoneToZoneInterpolation for the calculation of
       // pointDistanceToIntersection
       //- and faceDistanceToIntersection
       //zoneToZoneInterpolation* masterToSlaveInterpolatorPtr_;
       PatchToPatchInterpolation<
           PrimitivePatch<
               face, List, pointField
               >, PrimitivePatch<face, List, pointField>
           >* masterToSlaveInterpolatorPtr_;

    // Private Member Functions

        //- Disallow copy construct
        normalContactModel(const normalContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const normalContactModel&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return patch_.boundaryMesh().mesh();
        }


public:

    //- Runtime type information
    TypeName("normalContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            normalContactModel,
            dictionary,
            (
                const word name,
                const fvPatch& patch,
                const dictionary& dict,
                const label masterPatchID,
                const label slavePatchID,
                const label masterFaceZoneID,
                const label slaveFaceZoneID,
                const PrimitivePatch<
                face, List, pointField
                >& masterFaceZonePatch,
                const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
            ),
            (
                name,
                patch,
                dict,
                masterPatchID,
                slavePatchID,
                masterFaceZoneID,
                slaveFaceZoneID,
                masterFaceZonePatch,
                slaveFaceZonePatch
                )
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<normalContactModel> New
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID,
            const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
        );


    // Constructors

        //- Construct from dictionary
        normalContactModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID,
            const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
        );


    // Destructor

        virtual ~normalContactModel()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return reference to patch
        const fvPatch& patch() const
        {
            return patch_;
        }

        //- Correct contact model
        virtual void correct
        (
            const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
            const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
            const intersection::algorithm alg,
            const intersection::direction dir,
            word fieldName,
            const Switch orthotropic,
            const word nonLinear,
            vectorField& slaveFaceNormals,
            GGIInterpolation< PrimitivePatch< face, List, pointField >,
                              PrimitivePatch< face, List, pointField >
                              >* ggiInterpolatorPtr
            ) = 0;

        //- Return slave displacement
        virtual const vectorField& slaveDisp() const = 0;

        //- Return slave pressure
        virtual const vectorField& slavePressure() const = 0;

        //- Return slave valueFraction
        virtual const symmTensorField& slaveValueFrac() const = 0;

        //- Return master patch ID
        virtual label masterPatchID() const
        {
            return masterPatchID_;
        }

        //- Return master patch ID
        virtual label slavePatchID() const
        {
            return slavePatchID_;
        }

        //- Return master face zone ID
        virtual label masterFaceZoneID() const
        {
            return masterFaceZoneID_;
        }

        //- Return master face zone ID
        virtual label slaveFaceZoneID() const
        {
            return slaveFaceZoneID_;
        }

        //- Return point gap field
        virtual scalarField& slaveContactPointGap()
        {
            return slaveContactPointGap_;
        }

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
