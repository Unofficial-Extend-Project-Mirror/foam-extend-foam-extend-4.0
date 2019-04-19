/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::openFoamTableReader

Description
    Reads an interpolation table from a file - OpenFOAM-format

SourceFiles
    openFoamTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef openFoamTableReader_H
#define openFoamTableReader_H

#include "tableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class openFoamTableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class openFoamTableReader
:
    public tableReader<Type>
{

public:

    //- Runtime type information
    TypeName("openFoam");


    // Constructors

        //- Construct from dictionary
        openFoamTableReader(const dictionary &dict);

        //- Construct and return a copy
        virtual autoPtr<tableReader<Type> > clone() const
        {
            return autoPtr<tableReader<Type> >
            (
                new openFoamTableReader<Type>
                (
                    *this
                )
            );
        }


    //- Destructor
    virtual ~openFoamTableReader();


    // Member functions

        //- Read the table
        virtual void operator()(const fileName&, List<Tuple2<scalar, Type> > &);

        //- Read 2D table
        virtual void operator()
        (
            const fileName&,
            List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "openFoamTableReader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
