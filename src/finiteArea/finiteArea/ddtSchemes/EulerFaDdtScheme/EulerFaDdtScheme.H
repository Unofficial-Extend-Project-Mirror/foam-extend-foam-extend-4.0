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
    EulerDdtScheme

Description
    Basic first-order Euler implicit/explicit ddt using only the current and
    previous time-step values.

SourceFiles
    EulerDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerFaDdtScheme_H
#define EulerFaDdtScheme_H

#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class EulerFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerFaDdtScheme
:
    public fa::faDdtScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        EulerFaDdtScheme(const EulerFaDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const EulerFaDdtScheme&);


public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from mesh
        EulerFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        EulerFaDdtScheme(const faMesh& mesh, Istream& is)
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
#   include "EulerFaDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
