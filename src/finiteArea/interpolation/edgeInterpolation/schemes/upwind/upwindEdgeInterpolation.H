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
    upwindEdgeInterpolation

Description
    Upwind differencing scheme class.

SourceFiles
    upwindEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef upwindEdgeInterpolation_H
#define upwindEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private data

        const edgeScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const upwindEdgeInterpolation&);


public:

    //- Runtime type information
    TypeName("upwind");


    // Constructors

        //- Construct from faceFlux
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh().objectRegistry::lookupObject<edgeScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return pos(faceFlux_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
