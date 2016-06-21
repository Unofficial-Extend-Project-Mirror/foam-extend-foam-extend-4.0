/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::NoPostProcessing

Description
    Place holder for 'none' option

SourceFiles
    NoPostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef NoPostProcessing_H
#define NoPostProcessing_H

#include "PostProcessingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NoPostProcessing Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPostProcessing
:
    public PostProcessingModel<CloudType>
{
protected:

    // Protected member functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        NoPostProcessing(const dictionary&, CloudType&);


    //- Destructor
    virtual ~NoPostProcessing();


    // Member Functions

        // Evaluation

            //- Flag to indicate whether model activates post-processing model
            virtual bool active() const;

            //- Gather post-processing info on patch
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const label patchI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoPostProcessing.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
