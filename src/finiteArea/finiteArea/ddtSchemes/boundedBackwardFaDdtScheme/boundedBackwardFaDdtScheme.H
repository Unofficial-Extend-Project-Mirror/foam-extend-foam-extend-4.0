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
    boundedBackwardFaDdtScheme

Description
    Second-order backward-differencing implicit ddt using the current and
    two previous time-step values.

SourceFiles
    boundedBackwardFaDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef boundedBackwardFaDdtScheme_H
#define boundedBackwardFaDdtScheme_H

#include "faDdtScheme.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class boundedBackwardFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

class boundedBackwardFaDdtScheme
:
    public fa::faDdtScheme<scalar>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField& vf) const
        {
            if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
            {
                return GREAT;
            }
            else
            {
                return deltaT0_();
            }
        }


        //- Disallow default bitwise copy construct
        boundedBackwardFaDdtScheme(const boundedBackwardFaDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const boundedBackwardFaDdtScheme&);


public:

    //- Runtime type information
    TypeName("boundedBackward");


    // Constructors

        //- Construct from mesh
        boundedBackwardFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<scalar>(mesh)
        {}

        //- Construct from mesh and Istream
        boundedBackwardFaDdtScheme(const faMesh& mesh, Istream& is)
        :
            faDdtScheme<scalar>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return fa::faDdtScheme<scalar>::mesh();
        }

        tmp<areaScalarField> facDdt
        (
            const dimensionedScalar
        );

        tmp<areaScalarField> facDdt0
        (
            const dimensionedScalar
        );

        tmp<areaScalarField> facDdt
        (
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const areaScalarField&
        );

        tmp<edgeScalarField> facDdt0
        (
            const edgeScalarField&
        );

        tmp<areaScalarField> facDdt
        (
            const dimensionedScalar&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const dimensionedScalar&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt
        (
            const areaScalarField&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const areaScalarField&,
            const areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            const dimensionedScalar&,
            areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            const areaScalarField&,
            areaScalarField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
