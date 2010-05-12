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
    backwardFaDdtScheme

Description
    Second-order backward-differencing implicit ddt using the current and
    two previous time-step values.

SourceFiles
    backwardFaDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef backwardFaDdtScheme_H
#define backwardFaDdtScheme_H

#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class backwardFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class backwardFaDdtScheme
:
    public fa::faDdtScheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;

        //- Disallow default bitwise copy construct
        backwardFaDdtScheme(const backwardFaDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const backwardFaDdtScheme&);


public:

    //- Runtime type information
    TypeName("backward");


    // Constructors

        //- Construct from mesh
        backwardFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        backwardFaDdtScheme(const faMesh& mesh, Istream& is)
        :
            faDdtScheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return fa::faDdtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faePatchField, edgeMesh> > facDdt0
        (
            const GeometricField<Type, faePatchField, edgeMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh> > facDdt0
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type> > famDdt
        (
            GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type> > famDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type> > famDdt
        (
            const areaScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "backwardFaDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
