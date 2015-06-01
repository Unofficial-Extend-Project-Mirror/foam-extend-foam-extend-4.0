/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::anisotropicFilter

Description
    anisotropic filter

    @verbatim
    Kernel                 as filter          as Test filter with ratio 2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Box filter:            g = delta2/24  ->  g = delta2/6
    Spherical box filter:  g = delta2/64  ->  g = delta2/16
    Gaussian filter:       g = delta2/24  ->  g = delta2/6
    @endverbatim

SourceFiles
    anisotropicFilter.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicFilter_H
#define anisotropicFilter_H

#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class anisotropicFilter Declaration
\*---------------------------------------------------------------------------*/

class anisotropicFilter
:
    public LESfilter
{
    // Private data

        scalar widthCoeff_;
        volVectorField coeff_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        anisotropicFilter(const anisotropicFilter&);
        void operator=(const anisotropicFilter&);


public:

    //- Runtime type information
    TypeName("anisotropic");

    // Constructors

        //- Construct from components
        anisotropicFilter(const fvMesh& mesh, scalar widthCoeff);

        //- Construct from IOdictionary
        anisotropicFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~anisotropicFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()(const tmp<volScalarField>&) const;

        virtual tmp<volVectorField> operator()(const tmp<volVectorField>&) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()(const tmp<volTensorField>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
