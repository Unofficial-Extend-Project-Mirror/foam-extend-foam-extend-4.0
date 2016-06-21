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
    MixingPlaneInterpolation

Description
    The mixingPlaneInterplation class is implementing the interpolation of face
    data between two primitivePatches using mixing plane averaging.

    Some overview of this development can be found here:

    OpenFOAM TURBO TOOLS: FROM GENERAL PURPOSE CFD TO TURBOMACHINERY
    SIMULATIONS
    H. Jasak, M. Beaudoin, Proceedings of ASME-JSME-KSME Joint Fluids
    Engineering Conference 2011, AJK2011-FED, July 24-29, 2011, Hamamatsu,
    Shizuoka, JAPAN

    Steady-state capabilities for hydroturbines with OpenFOAM,
    M. Page, M. Beaudoin, A.-M. Giroux, Internationnal Journal of Fluid
    Machinery and Systems, 4(1):160–170, Jan-Mar 2011.

    Development of a General Grid Interface for Turbomachinery simulations with
    OpenFOAM, M. Beaudoin, H. Jasak, Open Source CFD International Conference,
    December 2008

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.


Giving credit where credit is due:
    1: Hakan Nilsson from the Chalmers University of Technology came up with
        the initial idea of using a combination of 2 GGI interfaces sharing a
        common single-face 360 degree ribbons patch in order to compute the
        circumferential average of fields. This implementation of the
        mixingPlane interpolation algorithm is an exploration of this simple
        but rather powerful idea.

    2: Maryse Page from Hydro-Quebec provided many test cases and many
       simulation runs for testing this interpolation algorithm. Testing is
       obviously an essential part of any new development.

    3: The authors also want to acknowledge the useful comments from many
       colleagues in the OpenFOAM Turbomachinery Special Interest Group.


Nomenclature
    ribbon = face object over which we are averaging

    direction = axis in which ribbon width is measured

    span = axis in which ribbon length is measured
           All ribbons have identical length in the chosen coordinate system

SourceFiles
    MixingPlaneInterpolation.C
    MixingPlaneProfile.C
    MixingPlaneInterpolationPatches.C
    MixingPlaneInterpolationAddressing.C
    MixingPlaneInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef MixingPlaneInterpolation_H
#define MixingPlaneInterpolation_H

#include "Pstream.H"
#include "className.H"
#include "labelList.H"
#include "scalarField.H"
#include "pointField.H"
#include "FieldFields.H"
#include "faceList.H"
#include "intersection.H"
#include "point2D.H"
#include "NamedEnum.H"
#include "coordinateSystem.H"
#include "boundBox.H"
#include "ggiInterpolation.H"
#include "cylindricalCS.H"
#include "primitivePatch.H"
#include "standAlonePatch.H"
#include <map>
#include <list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MixingPlaneInterpolationName Declaration
\*---------------------------------------------------------------------------*/

class MixingPlaneInterpolationName
{
public:

    // Public enumerations

        //- Define profile discretisation rule
        enum discretisation
        {
            MASTER_PATCH,
            SLAVE_PATCH,
            BOTH_PATCHES,
            UNIFORM,
            USER_DEFINED // User-defined through a file
        };

        //- Define the sweepAxis used to generate the mixing ribbon patch
        enum sweepAxis
        {
            // Cartesian coordinate system
            SWEEP_X,
            SWEEP_Y,
            SWEEP_Z,

            // Cylindrical coordinate system
            SWEEP_R,
            SWEEP_THETA,
            // SWEEP_Z is already defined for the Cartesian coordinate system
            SWEEP_UNKNOWN
        };

        //- Define the stackAxis used to generate the mixing ribbon patch
        enum stackAxis
        {
            // Cartesian coordinate system
            STACK_X,
            STACK_Y,
            STACK_Z,

            // Cylindrical coordinate system
            STACK_R,
            STACK_THETA,
            // STACK_Z is already defined for the Cartesian coordinate system
            STACK_UNKNOWN
        };

        //- Define type of mixing for field over patch
        enum mixingType
        {
            AREA_AVERAGING,
            FLUX_AVERAGING,
            UNIFORM_VALUE,
            UNIFORM_GRADIENT,
            ZERO_GRADIENT,
            MIXING_UNKNOWN
        };


    // Static data

        ClassName("MixingPlaneInterpolation");


       //- Discretisation names
       static const NamedEnum<discretisation, 5> discretisationNames_;

       //- Sweep axis names
       static const NamedEnum<sweepAxis, 6> sweepAxisNames_;

       //- Stack axis names
       static const NamedEnum<stackAxis, 6> stackAxisNames_;

       //- Mixing names
       static const NamedEnum<mixingType, 6> mixingTypeNames_;


    // Constructors

        //- Construct null
        MixingPlaneInterpolationName()
        {}
};


/*---------------------------------------------------------------------------*\
                    Class MixingPlaneInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class MasterPatch, class SlavePatch>
class MixingPlaneInterpolation
:
    public MixingPlaneInterpolationName
{
    // Private data types

    typedef std::map<Foam::scalar, std::list<Foam::point> > profileHistogram;


    // Private data

        //- Reference to the master patch
        const MasterPatch& masterPatch_;

        //- Reference to the slave patch
        const SlavePatch& slavePatch_;

        // We need to get some patches identification for debugging purposes
        // Since we are currently using primitive patches for the interpolator
        // we need to memorize this information when using the constructor

        //- Coordinate system for averaging
        const coordinateSystem& cs_;

        //- Type of mixing plane discretisation algorithm
        discretisation discretisationType_;

        //- Orientation of the mixing plane ribbon patch
        //  sweepAxis = axis in which ribbon length is measured
        //  stackAxis = axis in which ribbon width/height is measured
        sweepAxis sweepAxisType_;
        stackAxis stackAxisType_;

        //- Interpolation profile
        //  This list of points defines the profile generating the 'n'
        //  circular bands we are going to use for
        //  the circumferential averaging algorithm. For for 'n'
        //  circular bands, we need 'n + 1' points.
        mutable pointField interpolationProfile_;

        //- Slave-to-master transformation tensor.  Transforms slave data to
        //  master plane.  Size equals number of slave faces; zero length
        //  indicates no transform.  Size 1 indicates constant transform
        const tensorField forwardT_;

        //- Master-to-slave transformation tensor.  Transforms slave data to
        //  master plane.  Size equals number of master faces; zero length
        //  indicates no transform.  Size 1 indicates constant transform
        const tensorField reverseT_;

        //- Slave-to-master separation vector. Translation of slave data to
        //  master plane.  Size equals number of slave faces; zero length
        //  indicates no translation.
        const vectorField forwardSep_;


    // Demand-driven data

        //- Transformed master patch
        mutable standAlonePatch* transformedMasterPatchPtr_;

        //- Transformed shadow patch
        mutable standAlonePatch* transformedShadowPatchPtr_;

        //- Interpolating patch: strips of the mixing plane
        mutable standAlonePatch* mixingPlanePatchPtr_;


        //- Tensors for transforming the fields onto a single profile
        //  before averaging
        mutable tensorField* masterPatchToProfileTPtr_;
        mutable tensorField* masterProfileToPatchTPtr_;
        mutable tensorField* slavePatchToProfileTPtr_;
        mutable tensorField* slaveProfileToPatchTPtr_;


        // Master patch-to-profile interpolation

            //- Master addressing
            mutable labelListList* masterPatchToProfileAddrPtr_;
            mutable labelListList* masterProfileToPatchAddrPtr_;

            //- Master weights
            mutable scalarListList* masterPatchToProfileWeightsPtr_;
            mutable scalarListList* masterProfileToPatchWeightsPtr_;


        // Slave patch-to-profile interpolation

            //- Slave addressing
            mutable labelListList* slavePatchToProfileAddrPtr_;
            mutable labelListList* slaveProfileToPatchAddrPtr_;

            //- Slave weights
            mutable scalarListList* slavePatchToProfileWeightsPtr_;
            mutable scalarListList* slaveProfileToPatchWeightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MixingPlaneInterpolation
        (
            const MixingPlaneInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const MixingPlaneInterpolation&);


    // Helper functions for demand-driven data

        //- Return ribbon sweep direction
        direction sweepAxisSwitch() const;

        //- Return ribbon stack direction
        direction stackAxisSwitch() const;

        //- Construct profile histogram
        void updateProfileHistogram
        (
            profileHistogram& histo,
            const point& profileCoord,  // 3D point
            const direction dir,        // Sorting dimension 0: x, 1: y, 2: z
            const scalar halfSizeBin    // half size of min histogram bin width
        ) const;

        //- Create an interpolation profile from histograms
        //  Note: histograms are adjusted to achieve full overlap
        //  HJ, 6/Nov/2009
        tmp<pointField> computeProfileFromHistograms
        (
            const profileHistogram& masterHisto,
            const profileHistogram& slaveHisto,
            const scalar halfSizeBin   // half size of min histogram bin width
        ) const;

        //- Remove interpolationProfile_ points located outside of either
        //  master/slave patch boundingBox,
        //  with the exception of the first and last profile points
        void removeNonOverlappedProfilePoints
        (
            boundBox& masterBB,
            boundBox& slaveBB
        ) const;

        //- Find and modify patch faces straddling the histogram span
        //  eg. -180,+180 angle cylindrical coordinate system
        void correctStraddlingFaces
        (
            faceList& cylCoordFaces,
            pointField& cylCoordFacesPoint
        ) const;


        //- Calculate transformed patches
        void calcTransformedPatches() const;

        //- Calculate mixing patch
        void calcMixingPlanePatch() const;


        //- Calculate addressing and weights
        void calcAddressing() const;

        //- Calculate Cartesian to cylindrical tranformation tensor field for
        //  master and slave patches
        void calcTransforms() const;

        //- Create an interpolation profile from histograms
        tmp<pointField> calcProfile() const;


        //- Clear transformed patches
        void clearTransfomedPatches();

        //- Calculate mixing patch
        void clearMixingPlanePatch();

        //- Clear addressing
        void clearAddressing();

        //- Clear transforms
        void clearTransforms();

        //- Clear all geometry and addressing
        void clearOut();


    // Interpolation data access

        //- Distribute to profile given addressing and weights
        template<class Type>
        static void toProfile
        (
            const Field<Type>& srcF,
            const labelListList& srcAddr,
            const scalarListList& srcWeights,
            Field<Type>& profileBandValues
        );

        //- Collect from profile given addressing and weights
        template<class Type>
        static void fromProfile
        (
            const Field<Type>& profileBandValues,
            const labelListList& dstAddr,
            const scalarListList& dstWeights,
            Field<Type>& dstResultF
        );

        //- Collect from profile given addressing and weights
        //  for masked faces only
        template<class Type>
        static void maskedFromProfile
        (
            const Field<Type>& profileBandValues,
            const labelListList& dstAddr,
            const scalarListList& dstWeights,
            Field<Type>& dstResultF,
            const labelList& mask
        );

        //- Execute transform for masked faces only
        template<class Type>
        static void maskedTransform
        (
            Field<Type>& transField,
            const tensorField& t,
            const Field<Type>& inField,
            const labelList& mask
        );

        //- Interpolate given source and target, addressing and weights
        template<class Type>
        void interpolate
        (
            const Field<Type>& srcF,
            const labelListList& srcAddr,
            const scalarListList& srcWeights,
            const labelListList& dstAddr,
            const scalarListList& dstWeights,
            Field<Type>& dstResultF
        ) const;


        //- Is a transform required?
        inline bool doTransform() const
        {
            return false;
        }

        //- Is a separation required?
        inline bool doSeparation() const
        {
            return false;
        }


public:

    // Constructors

        //- Construct from components
        MixingPlaneInterpolation
        (
            const MasterPatch& masterPatch,
            const SlavePatch& slavePatch,
            const coordinateSystem& cs,
            const discretisation& discretisationType,
            const sweepAxis& sweepAxisType,
            const stackAxis& stackAxisType,
            const pointField& interpolationProfile
        );


    // Destructor

        ~MixingPlaneInterpolation();


    // Member Functions

        // Access

            //- Return master patch
            const MasterPatch& masterPatch() const
            {
                return masterPatch_;
            }

            //- Return slave patch
            const SlavePatch& slavePatch() const
            {
                return slavePatch_;
            }

            //- Return number of profile bands
            label nProfileBands() const
            {
                return interpolationProfile().size() - 1;
            }

            //- Return interpolation profile
            const pointField& interpolationProfile() const;


            //- Return mixing plane patch
            const standAlonePatch& mixingPlanePatch() const;

            //- Return transformed patch
            const standAlonePatch& transformedMasterPatch() const;

            //- Return transformed patch
            const standAlonePatch& transformedShadowPatch() const;


            //- Return reference to master addressing
            const labelListList& masterPatchToProfileAddr() const;
            const labelListList& masterProfileToPatchAddr() const;

            //- Return reference to master weights
            const scalarListList& masterPatchToProfileWeights() const;
            const scalarListList& masterProfileToPatchWeights() const;

            //- Return reference to slave addressing
            const labelListList& slavePatchToProfileAddr() const;
            const labelListList& slaveProfileToPatchAddr() const;

            //- Return reference to slave weights
            const scalarListList& slavePatchToProfileWeights() const;
            const scalarListList& slaveProfileToPatchWeights() const;

            //- Return reference to masterPatchToProfile tensor field
            const tensorField& masterPatchToProfileT() const;

            //- Return reference to masterProfileToPatchT tensor field
            const tensorField& masterProfileToPatchT() const;

            //- Return reference to slavePatchToProfileT tensor field
            const tensorField& slavePatchToProfileT() const;

            //- Return reference to slaveProfileToPatchT tensor field
            const tensorField& slaveProfileToPatchT() const;


        // Interpolation functions

            //- Interpolate from master to slave
            template<class Type>
            tmp<Field<Type> > masterToSlave
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > masterToSlave
            (
                const tmp<Field<Type> >& tpf
            ) const;


            //- Interpolate from slave to master
            template<class Type>
            tmp<Field<Type> > slaveToMaster
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > slaveToMaster
            (
                const tmp<Field<Type> >& tpf
            ) const;


            //- Interpolate from master to profile
            template<class Type>
            tmp<Field<Type> > masterToProfile
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > masterToProfile
            (
                const tmp<Field<Type> >& tpf
            ) const;

            //- Interpolate from slave to profile
            template<class Type>
            tmp<Field<Type> > slaveToProfile
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > slaveToProfile
            (
                const tmp<Field<Type> >& tpf
            ) const;


            //- Interpolate from profile to master
            template<class Type>
            tmp<Field<Type> > profileToMaster
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > profileToMaster
            (
                const tmp<Field<Type> >& tpf
            ) const;

            template<class Type>
            void maskedProfileToMaster
            (
                const Field<Type>& profileFF,
                Field<Type>& result,
                const labelList& mask
            ) const;

            //- Interpolate from profile to slave
            template<class Type>
            tmp<Field<Type> > profileToSlave
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > profileToSlave
            (
                const tmp<Field<Type> >& tpf
            ) const;


            template<class Type>
            void maskedProfileToSlave
            (
                const Field<Type>& profileFF,
                Field<Type>& result,
                const labelList& mask
            ) const;

            //- Interpolate from master to master
            //- (self circumferential averaging)
            template<class Type>
            tmp<Field<Type> > masterToMaster
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > masterToMaster
            (
                const tmp<Field<Type> >& tpf
            ) const;

            //- Interpolate from slave to slave
            //- (self circumferential averaging)
            template<class Type>
            tmp<Field<Type> > slaveToSlave
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > slaveToSlave
            (
                const tmp<Field<Type> >& tpf
            ) const;


        // Edit

            //- Correct weighting factors for moving mesh.
            bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "MixingPlaneInterpolationTemplate.C"
#   include "MixingPlaneProfile.C"
#   include "MixingPlaneInterpolationPatches.C"
#   include "MixingPlaneInterpolationAddressing.C"
#   include "MixingPlaneInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
