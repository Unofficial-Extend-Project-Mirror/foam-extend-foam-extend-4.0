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
    frictionContactModel

Description
    Model for the the contact friction
    The actual friction law (e.g. Coulomb) is specified by the frictionLaw
    in the frictionLaw

SourceFiles
    frictionContactModel.C
    newFrictionContactModel.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef frictionContactModel_H
#define frictionContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class frictionContactModel Declaration
\*---------------------------------------------------------------------------*/

class frictionContactModel
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

        // for visualistion of sticking and slipping faces
        // face not in contact is set to 0.0
        // slipping face is set to 1.0
        // sticking face is set to 2.0
        scalarField stickSlipFaces_;


    // Private Member Functions

        //- Disallow copy construct
        frictionContactModel(const frictionContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const frictionContactModel&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
          return patch_.boundaryMesh().mesh();
        }


public:

    //- Runtime type information
    TypeName("frictionContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionContactModel,
            dictionary,
            (
                const word name,
                const fvPatch& patch,
                const dictionary& dict,
                const label masterPatchID,
                const label slavePatchID,
                const label masterFaceZoneID,
                const label slaveFaceZoneID
            ),
            (
                name,
                patch,
                dict,
                masterPatchID,
                slavePatchID,
                masterFaceZoneID,
                slaveFaceZoneID
            )
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<frictionContactModel> New
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const label masterFaceZoneID,
            const label slaveFaceZoneID
        );


    // Constructors

        //- Construct from dictionary
        frictionContactModel
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

        virtual ~frictionContactModel()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Correct contatc model
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
        ) = 0;

        //- Return slave friction displacement
        virtual const vectorField& slaveDisp() const = 0;

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const = 0;

        //- Return slave friction value fraction
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

        //- Return stick slip faces field
        // virtual volScalarField& stickSlipFaces()
        virtual scalarField& stickSlipFaces()
        {
            return stickSlipFaces_;
        }

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
