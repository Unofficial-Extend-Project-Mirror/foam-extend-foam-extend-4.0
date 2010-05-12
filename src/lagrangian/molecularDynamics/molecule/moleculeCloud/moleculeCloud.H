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
    Foam::moleculeCloud

Description

SourceFiles
    moleculeCloudApplyConstraintsAndThermostats.C
    moleculeCloudBuildCellInteractionLists.C
    moleculeCloudBuildCellOccupancy.C
    moleculeCloudBuildCellReferralLists.C
    moleculeCloud.C
    moleculeCloudCalculateAndAccumulateProperties.C
    moleculeCloudCalculateExternalForce.C
    moleculeCloudCalculateForce.C
    moleculeCloudCalculatePairForce.C
    moleculeCloudCalculateTetherForce.C
    moleculeCloudIntegrateEquationsOfMotion.C
    moleculeCloudRemoveHighEnergyOverlaps.C


\*---------------------------------------------------------------------------*/

#ifndef moleculeCloud_H
#define moleculeCloud_H

#include "Cloud.H"
#include "molecule.H"
#include "IOdictionary.H"
#include "vector2D.H"

#include "pairPotentialList.H"
#include "tetherPotentialList.H"

#include "receivingReferralList.H"
#include "sendingReferralList.H"
#include "referredCellList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class moleculeCloud Declaration
\*---------------------------------------------------------------------------*/

class moleculeCloud
:
    public Cloud<molecule>
{

public:

    enum integrationMethods
    {
        imVerletLeapfrog,
        imPredictorCorrector
    };


private:

    // Private data

        const polyMesh& mesh_;

            // MD solution parameters

            static const NamedEnum<integrationMethods, 2>
                integrationMethodNames_;

            integrationMethods integrationMethod_;

            scalar potentialEnergyLimit_;

            labelList removalOrder_;

            labelListList directInteractionList_;

            referredCellList referredInteractionList_;

            labelList realCellsWithinRCutMaxOfAnyReferringPatch_;

            labelList realFacesWithinRCutMaxOfAnyReferringPatch_;

            labelList realEdgesWithinRCutMaxOfAnyReferringPatch_;

            labelList realPointsWithinRCutMaxOfAnyReferringPatch_;

            List<sendingReferralList> cellSendingReferralLists_;

            List<receivingReferralList> cellReceivingReferralLists_;

            pairPotentialList pairPotentials_;

            tetherPotentialList tetherPotentials_;

            vector gravity_;

            List< DynamicList<molecule*> > cellOccupancy_;


        // Private Member Functions

            //- Disallow default bitwise copy construct
            moleculeCloud(const moleculeCloud&);

            //- Disallow default bitwise assignment
            void operator=(const moleculeCloud&);


public:

    // Static data members

        //- Tolerance for checking that faces on a patch segment
        static scalar transTol;

        static scalar kb;


    // Constructors

        //- Construct given mesh
        moleculeCloud(const polyMesh&);

        //- Construct given polyMesh and fields of position, cell, mass,
        //- id, U ands A.  Intended for use by the molConfig utility
        moleculeCloud
        (
            const polyMesh& mesh,
            label nMol,
            const labelField& id,
            const scalarField& mass,
            const vectorField& positions,
            const labelField& cells,
            const vectorField& U,
            const vectorField& A,
            const labelField& tethered,
            const vectorField& tetherPositions
        );


    // Member Functions

        // Access

            inline const polyMesh& mesh() const;

        // MD solution parameters

            inline const integrationMethods& integrationMethod() const;

            inline scalar potentialEnergyLimit() const;

            inline const labelList& removalOrder() const;

            inline label nPairPotentials() const;

            inline const labelListList& directInteractionList() const;

            inline const referredCellList& referredInteractionList() const;

            inline const labelList&
                realCellsWithinRCutMaxOfAnyReferringPatch() const;

            inline const labelList&
                realFacesWithinRCutMaxOfAnyReferringPatch() const;

            inline const labelList&
                realEdgesWithinRCutMaxOfAnyReferringPatch() const;

            inline const labelList&
                realPointsWithinRCutMaxOfAnyReferringPatch() const;

            inline const List<sendingReferralList>&
                cellSendingReferralLists() const;

            inline const List<receivingReferralList>&
                cellReceivingReferralLists() const;

            inline label nInteractingProcs() const;

            inline const pairPotentialList& pairPotentials() const;

            inline const tetherPotentialList& tetherPotentials() const;

            inline const vector& gravity() const;

            inline const List< DynamicList<molecule*> >& cellOccupancy() const;

            void buildCellInteractionLists();

            //- Build referralLists which define how to send information
            // to referredCells to source cells
            void buildCellReferralLists();

            bool testPointFaceDistance
            (
                const label p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const label p,
                const referredCell& refCell
            ) const;

            bool testPointFaceDistance
            (
                const vectorList& pointsToTest,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const labelList& faceToTest,
                const vectorList& points,
                const vector& faceC,
                const vector& faceA
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const edge& eJ
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const vector& eJs,
                const vector& eJe
            ) const;

            const labelList realCellsInRangeOfSegment
            (
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;

            const labelList referredCellsInRangeOfSegment
            (
                const List<referredCell>& referredInteractionList,
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;

            //- Determine which molecules are in which cells
            void buildCellOccupancy();

            //- Integrate the equations of motion using algorithm selected at
            //  runtime from a dictionary. This will also call the function
            //  to calculate the intermolecular forces (calculatePairForce()).

            void integrateEquationsOfMotion();

            void applyConstraintsAndThermostats
            (
                const scalar targetTemperature,
                const scalar measuredTemperature
            );

            void calculateForce();

            void calculatePairForce();

            void calculateTetherForce();

            void calculateExternalForce();

            void removeHighEnergyOverlaps();


    // Member Operators

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
