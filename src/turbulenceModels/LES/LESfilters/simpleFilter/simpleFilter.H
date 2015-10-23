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
    Foam::simpleFilter

Description
    Simple top-hat filter used in dynamic LES models.

    Implemented as a surface integral of the face interpolate of the field.

SourceFiles
    simpleFilter.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFilter_H
#define simpleFilter_H

#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleFilter Declaration
\*---------------------------------------------------------------------------*/

class simpleFilter
:
    public LESfilter

{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        simpleFilter(const simpleFilter&);
        void operator=(const simpleFilter&);


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct from components
        simpleFilter(const fvMesh& mesh);

        //- Construct from IOdictionary
        simpleFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~simpleFilter()
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
