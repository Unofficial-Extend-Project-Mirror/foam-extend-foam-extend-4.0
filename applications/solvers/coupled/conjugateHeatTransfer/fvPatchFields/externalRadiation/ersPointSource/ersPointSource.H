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
    Foam::ersPointSource

Description
    helper class to sum up external radiation sources


SourceFiles
    ersPointSource.C

\*---------------------------------------------------------------------------*/

#ifndef ersPointSource_H
#define ersPointSource_H

#include "ersConstantFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wallHeatTransferFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ersPointSource
:
    public ersConstantFlux
{
    // Private data

        //- Strenght of the point source
        scalar qmax_;

        //- radiation absorptivity
        scalar alpha_;

        //- Direction of the ersPointSource
        vector direction_;

public:

    //- Runtime type information
    TypeName("pointSource");


    // Constructors

        //- Construct from components
        explicit ersPointSource
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& p
        );


        //- Construct and return a clone
        virtual autoPtr<externalRadiationSource> clone() const
        {
            return autoPtr<externalRadiationSource>
            (
                new ersPointSource(*this)
            );
        }

    // Destructor

        virtual ~ersPointSource()
        {};


    // Member functions

        // Access

        // Mapping functions

        // Evaluation functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
