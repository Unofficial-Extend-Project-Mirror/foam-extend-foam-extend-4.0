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
    Foam::sampledSurface

Description
    An abstract class for surfaces with sampling.

SourceFiles
    sampledSurface.C
    sampledSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurface_H
#define sampledSurface_H

#include "pointField.H"
#include "word.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "polyMesh.H"
#include "coordinateSystems.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sampledSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledSurface
{
    // Private data

        //- Name of sample surface
        word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Triangulate faces or keep faces as it
        bool triangulate_;

        //- Do we intend to interpolate the information?
        bool interpolate_;


    // Demand-driven data

        //- Face area vectors
        mutable vectorField* SfPtr_;

        //- Mag face area vectors
        mutable scalarField* magSfPtr_;

        //- Face centres
        mutable vectorField* CfPtr_;

        //- Total surface area
        mutable scalar area_;

    // Make geometric data

        //- Make Sf
        void makeSf() const;

        //- Make magSf
        void makeMagSf() const;

        //- Make Cf
        void makeCf() const;

    // Service methods

        //- Check field size matches surface size
        template<class Type>
        bool checkFieldSize(const Field<Type>&) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        void project
        (
            Field<ReturnType>&,
            const Field<Type>&
        ) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        void project
        (
            Field<ReturnType>&,
            const tmp<Field<Type> >&
        ) const;

        //- Project field onto surface
        template<class ReturnType, class Type>
        tmp<Field<ReturnType> > project(const tmp<Field<Type> >&) const;

protected:

    // Protected Member functions

        virtual void clearGeom() const;

        //- Non-const access to triangulation
        bool& triangulate()
        {
            return triangulate_;
        }

public:

    //- Runtime type information
      TypeName("sampledSurface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sampledSurface,
            word,
            (
                const word& name,
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


        //- Class used for the PtrLists read-construction
        //  Has the ability to rewrite coordinate systems as required
        class iNew
        :
            public coordinateSystems
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                coordinateSystems(mesh),
                mesh_(mesh)
            {}

            iNew
            (
                const polyMesh& mesh,
                const coordinateSystems& cs
            )
            :
                coordinateSystems(cs),
                mesh_(mesh)
            {}

            autoPtr<sampledSurface> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                rewriteDict(dict, true);

                return sampledSurface::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from name, mesh
        sampledSurface
        (
            const word& name,
            const polyMesh& mesh,
            const bool triangulate = true
        );

        //- Construct from dictionary
        sampledSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Clone
        autoPtr<sampledSurface> clone() const
        {
            notImplemented("autoPtr<sampledSurface> clone() const");
            return autoPtr<sampledSurface>(NULL);
        }


    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<sampledSurface> New
        (
            const word& sampleType,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~sampledSurface();


    // Member Functions

      // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Name of surface
        const word& name() const
        {
            return name_;
        }

        //- interpolation requested for surface
        bool interpolate() const
        {
            return interpolate_;
        }

        //- triangulation requested for surface
        bool triangulate() const
        {
            return triangulate_;
        }

        //- Points of surface
        virtual const pointField& points() const = 0;

        //- Faces of surface
        virtual const faceList& faces() const = 0;

        //- Correct for mesh movement and/or field changes
        virtual void correct(const bool meshChanged) = 0;

        //- Return face area vectors
        virtual const vectorField& Sf() const;

        //- Return face area magnitudes
        virtual const scalarField& magSf() const;

        //- Return face centres as vectorField
        virtual const vectorField& Cf() const;

        //- The total surface area
        scalar area() const;

        //- Integration of a field across the surface
        template<class Type>
        Type integrate(const Field<Type>&) const;

        //- Integration of a field across the surface
        template<class Type>
        Type integrate(const tmp<Field<Type> >&) const;

        //- Area-averaged value of a field across the surface
        template<class Type>
        Type average(const Field<Type>&) const;

        //- Area-averaged value of a field across the surface
        template<class Type>
        Type average(const tmp<Field<Type> >&) const;

        //- Project field onto surface
        tmp<Field<scalar> > project(const Field<scalar>&) const;

        //- Project field onto surface
        tmp<Field<scalar> > project(const Field<vector>&) const;

        //- Project field onto surface
        tmp<Field<vector> > project(const Field<sphericalTensor>&) const;

        //- Project field onto surface
        tmp<Field<vector> > project(const Field<symmTensor>&) const;

        //- Project field onto surface
        tmp<Field<vector> > project(const Field<tensor>&) const;

        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const = 0;

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const = 0;


        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const = 0;


        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const = 0;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const = 0;


      // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }


      // Write

        //- Write
        virtual void print(Ostream&) const;


       // IOstream operators

        friend Ostream& operator<<(Ostream&, const sampledSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sampledSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
