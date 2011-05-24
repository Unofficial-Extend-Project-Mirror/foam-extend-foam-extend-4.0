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
    DeferredCorrectionLimitedScheme

Description

    Class to create NVD/TVD limited weighting-factors.
    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled patches
    (e.g. processor-processor patches). The weight function is supplied the
    central-differencing weighting factor, the face-flux, the cell and face
    gradients (from which the normalised variable distribution may be created)
    and the cell centre distance.

    This code organisation is both neat and efficient, allowing for convenient
    implementation of new schemes to run on parallelised cases.

    The scheme is implemeted in deferred correction mode, with upwind weighting
    factors being implicit.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    DeferredCorrectionLimitedScheme.C

\*---------------------------------------------------------------------------*/

#ifndef DeferredCorrectionLimitedScheme_H
#define DeferredCorrectionLimitedScheme_H

#include "limitedSurfaceInterpolationScheme.H"
#include "LimitFuncs.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class DeferredCorrectionLimitedScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Limiter, template<class> class LimitFunc>
class DeferredCorrectionLimitedScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public Limiter
{
    // Private data

        //- Upwind scheme for implicit part
        upwind<Type> upwindScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        DeferredCorrectionLimitedScheme(const DeferredCorrectionLimitedScheme&);

        //- Disallow default bitwise assignment
        void operator=(const DeferredCorrectionLimitedScheme&);


public:

    //- Runtime type information
    TypeName("DeferredCorrectionLimitedScheme");

    typedef Limiter LimiterType;

    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        DeferredCorrectionLimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const Limiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(weight),
            upwindScheme_(mesh, faceFlux)
        {}

        //- Construct from mesh and Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        DeferredCorrectionLimitedScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is),
            Limiter(is),
            upwindScheme_(mesh, this->faceFlux_)
        {}

        //- Construct from mesh, faceFlux and Istream
        DeferredCorrectionLimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(is),
            upwindScheme_(mesh, faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return upwindScheme_.weights();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the constructor functions to the hash tables

#define makeDeferredSurfaceInterpolationTypeScheme(SS, LIMITER, NVDTVD, LIMFUNC, TYPE) \
                                                                              \
typedef DeferredCorrectionLimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>             \
    DeferredCorrectionLimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_;       \
defineTemplateTypeNameAndDebugWithName                                        \
    (DeferredCorrectionLimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);               \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<DeferredCorrectionLimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >\
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<DeferredCorrectionLimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >\
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                      \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable            \
<DeferredCorrectionLimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >\
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                   \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable        \
<DeferredCorrectionLimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >\
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


#define makeDeferredSurfaceInterpolationScheme(SS, LIMITER)                    \
                                                                              \
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,scalar)    \
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,vector)    \
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,sphericalTensor) \
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,symmTensor)\
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,tensor)


#define makeDeferredVSurfaceInterpolationScheme(SS, LIMITER)                   \
makeDeferredSurfaceInterpolationTypeScheme(SS,LIMITER,NVDVTVDV,null,vector)


#define makeLDeferredSurfaceInterpolationTypeScheme(SS, LLIMITER, LIMITER, NVDTVD, LIMFUNC, TYPE) \
                                                                              \
typedef DeferredCorrectionLimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC>  \
    DeferredCorrectionLimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;               \
defineTemplateTypeNameAndDebugWithName                                        \
    (DeferredCorrectionLimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0); \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<DeferredCorrectionLimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> > \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<DeferredCorrectionLimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> > \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                      \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable            \
<DeferredCorrectionLimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> > \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                   \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable        \
<DeferredCorrectionLimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> > \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DeferredCorrectionLimitedScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
