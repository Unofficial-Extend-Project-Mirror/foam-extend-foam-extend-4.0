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
    Foam::exponential

Description
    exponential pdf

SourceFiles
    exponentialI.H
    exponential.C
    exponentialIO.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public pdf
{
    // Private data

        dictionary pdfDict_;

        //- min and max values of the distribution
        scalar minValue_;
        scalar maxValue_;

        List<scalar> lambda_;
        List<scalar> ls_;

        scalar range_;

public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        //- Construct from components
        exponential
        (
            const dictionary& dict,
            Random& rndGen
        );


    //- Destructor
    virtual ~exponential();


    // Member Functions

         //- Sample the pdf
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "exponentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
