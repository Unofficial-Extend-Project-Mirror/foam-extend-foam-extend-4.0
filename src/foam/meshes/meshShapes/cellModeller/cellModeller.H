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
    Foam::cellModeller

Description
    A static collection of cell models, and a means of looking them up.

SourceFiles
    cellModeller.C
    cellModellerIO.C
    globalCellModeller.C

\*---------------------------------------------------------------------------*/

#ifndef cellModeller_H
#define cellModeller_H

#include "cellModel.H"
#include "PtrList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellModeller Declaration
\*---------------------------------------------------------------------------*/

class cellModeller
{
    // Private data

        //- PtrList of models
        static PtrList<cellModel> models_;

        //- List of model pointers
        static List<cellModel*> modelPtrs_;

        //- HashTable of model pointers
        static HashTable<const cellModel*> modelDictionary_;


public:

    // Constructors

        //- Construct from central "cellModels" file
        cellModeller();

    // Destructor

        ~cellModeller();


    // Member functions

        //- Look up a model by name and return a pointer to the model or NULL
        static const cellModel* lookup(const word&);

        //- Look up a model by index and return a pointer to the model or NULL
        static const cellModel* lookup(const label i)
        {
            return modelPtrs_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
