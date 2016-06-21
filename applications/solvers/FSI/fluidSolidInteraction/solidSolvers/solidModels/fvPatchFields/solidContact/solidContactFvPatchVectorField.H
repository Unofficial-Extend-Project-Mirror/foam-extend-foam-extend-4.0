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
    solidContactFvPatchVectorField

Description
    contact boundary conditions for finite volume solidMechanics solvers.
    There is a master and a slave patch. The master creates the contactModel
    and the slave gets the pointer to it.
    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

SourceFiles
    solidContactFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef solidContactFvPatchVectorField_H
#define solidContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "GGIInterpolationTemplate.H"
#include "Switch.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidContactFvPatchVectorField
:
        public directionMixedFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- is this patch the master
        const Switch master_;

        //- Shadow patch name
        word shadowPatchName_;

        //- Shadow patch index
        label shadowPatchIndex_;

        //- Zone index
        mutable label* zoneIndexPtr_;

        //- Shadow zone index
        mutable label* shadowZoneIndexPtr_;

        // is the contact boundary active
        Switch contactActive_;

        // is the contact boundary active
        Switch rigidMaster_;

        // Copy of contact dictionary
        const dictionary dict_;

        //- Normal contact model pointer
        mutable normalContactModel* normalModelPtr_;

        //- Friction contact model pointer
        mutable frictionContactModel* frictionModelPtr_;

        // Mesha.allPoints() where the contact zone points are kept in the
        // deformed configuration
        mutable vectorField* allPointsDeformedPtr_;

        // Zone points at end of previous time step for incremental solver,
        // zone points from initial configuration for total displacement solver
        mutable vectorField* oldZonePointsPtr_;

        // Shadow zone points at end of previous time step for incremental
        // solver,
        // shadow zone points from initial configuration for total displacement
        // solver
        mutable vectorField* oldShadowZonePointsPtr_;

        // Master face zone kept in deformed configuration
        mutable primitiveFacePatch* zonePtr_;

        // Slave face zone kept in deformed configuration
        mutable primitiveFacePatch* shadowZonePtr_;

        //- Zone-to-zone interpolation
        //PrimitivePatch<face, List, const pointField&>
        mutable ggiZoneInterpolation* zoneToZonePtr_;

        //- patch to patch algorithm
        const intersection::algorithm alg_;

        //- patch to patch direction
        const intersection::direction dir_;

       // current time index
       label curTimeIndex_;

       // current momentum iterator
       label iCorr_;

       // freqency the contact is corrected
       label correctionFreq_;

       // if the solver uses an orthotropic approach
       Switch orthotropic_;

       // Force correction of contact on next call
       bool forceCorrection_;

       // Force correction of contact on next call
       // just friction part
       bool forceCorrectionFriction_;

       // Slave face normals
       // currently calculated by normalContactModel
       mutable vectorField* slaveFaceNormalsPtr_;

        //- Global face zones in the mesh
        mutable labelList* globalFaceZonesPtr_;

       // Point map from global to local for global face zones
       mutable IOList<labelList>* globalToLocalFaceZonePointMapPtr_;

       // Contact residual norm
       scalar contactResidual_;

       //- Type of geometrical non-linearity
       //- This is used to decide how to move the contact face zones
       nonLinearGeometry::nonLinearType nonLinear_;

    // Private Member Functions

            //- Calculate the zone index
            void calcZoneIndex() const;

            //- Calculate the shadow zone index
            void calcShadowZoneIndex() const;

            //- Create the normal contact model
            void calcNormalModel() const;

            //- Create the friction contact model
            void calcFrictionModel() const;

            //- Calculate deformed all points
            void calcAllPointsDeformed() const;

            //- Update deformed all points to deformed position
            void updateAllPointsDeformed();

            //- Calculate zone old points
            void calcOldZonePoints() const;

            //- Calculate shadow zone old points
            void calcOldShadowZonePoints() const;

            //- Calculate zone
            void calcZone() const;

            //- Calculate shadow zone
            void calcShadowZone() const;

            //- Calculate interpolation
            void calcZoneToZone() const;

            //- Calculate slaveFaceNormals
            void calcSlaveFaceNormals() const;

            //- Calculate global face zones
            void calcGlobalFaceZones() const;

            //- Calculate global to local face zone point map
            void calcGlobalToLocalFaceZonePointMap() const;

            // Contact dictionary
            const dictionary& dict() const
            {
                return dict_;
            };

            //- Get zone field from patch field
            template<class Type>
            tmp<Field<Type> > zoneField
            (
                const label zoneIndex,
                const label patchIndex,
                const Field<Type>& patchField
            ) const;

            //- Get patch field from zone field
            template<class Type>
            tmp<Field<Type> > patchField
            (
                const label patchIndex,
                const label zoneIndex,
                const Field<Type>& zoneField
            ) const;

            //- Interpolate zone face values to points
            template<class Type>
            tmp<Field<Type> > zoneFaceToPointInterpolate
            (
                const label zoneIndex,
                const Field<Type>& zoneField
            ) const;

            //- Interpolate zone point values to faces
            template<class Type>
            tmp<Field<Type> > zonePointToFaceInterpolate
            (
                const label zoneIndex,
                const Field<Type>& zonePointField
            ) const;

public:

    //- Runtime type information
    TypeName("solidContact");


    // Constructors

        //- Construct from patch and internal field
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidContactFvPatchVectorField onto
        //  a new patch
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidContactFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidContactFvPatchVectorField();


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return if this patch is the master
            Switch master() const
            {
                return master_;
            }

            //- Return the name of the shadow patch
            const word& shadowPatchName() const
            {
                return shadowPatchName_;
            }

            //- Return the name of the shadow patch
            label shadowPatchIndex() const
            {
                return shadowPatchIndex_;
            }

            //- Return index of the zone
            label zoneIndex() const;

            //- Return index of the shadow zone
            label shadowZoneIndex() const;

            //- Return reference to the normal contact model
            normalContactModel& normalModel();

            //- Return const reference to the normal contact model
            const normalContactModel& normalModel() const;

            //- Return reference to the friction contact model
            frictionContactModel& frictionModel();

            //- Return const reference to the friction contact model
            const frictionContactModel& frictionModel() const;

            //- Return reference to all points with zone points in deformed
            // position
            vectorField& allPointsDeformed();

            //- Return reference to all points with zone points in deformed
            // position
            const vectorField& allPointsDeformed() const;

            //- Return reference to old zone points
            const vectorField& oldZonePoints() const;

            //- Return reference to old shadow zone points
            const vectorField& oldShadowZonePoints() const;

            //- Return const reference to deformed master zone
            const primitiveFacePatch& zone() const;

            //- Return reference to deformed master zone
            primitiveFacePatch& zone();

            //- Return const reference to deformed master zone
            const primitiveFacePatch& shadowZone() const;

            //- Return reference to deformed master zone
            primitiveFacePatch& shadowZone();

            //- Return const reference to patch-to-patch interpolation
            const ggiZoneInterpolation& zoneToZone() const;

            //- Return reference to patch-to-patch interpolation
            ggiZoneInterpolation& zoneToZone();

            //- Return reference to slave face normals
            vectorField& slaveFaceNormals();

            //- Return const reference to global face zones
            const labelList& globalFaceZones() const;

            //- Return const reference to global face zone point map
            const IOList<labelList>& globalToLocalFaceZonePointMap() const;

            //- Return correction frequency
            label correctionFreq() const;

            //  Return contactActive
            virtual Switch contactActive() const
            {
                return contactActive_;
            }


        // Mapping functions

            //- Clear out demand driven data
            virtual void clearOut();

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;


        // Member functions

            //- Force contact correction
            virtual void forceCorrection()
            {
                forceCorrection_ = true;
            }

            //- Force contact correction
            virtual void forceCorrectionFriction()
            {
                forceCorrectionFriction_ = true;
            }

            //- Dissipated energy due to friction
            virtual tmp<scalarField> Qc() const;

            //- Return contact residual
            virtual scalar contactResidual()
            {
                return contactResidual_;
            }

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "solidContactFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
