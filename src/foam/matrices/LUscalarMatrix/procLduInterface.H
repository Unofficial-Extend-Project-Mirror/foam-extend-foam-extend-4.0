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
    Foam::procLduInterface

Description
    Foam::procLduInterface

SourceFiles
    procLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef procLduInterface_H
#define procLduInterface_H

#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduInterfaceField;

/*---------------------------------------------------------------------------*\
                           Class procLduInterface Declaration
\*---------------------------------------------------------------------------*/

class procLduInterface
{
    // Private data

        labelList faceCells_;
        scalarField coeffs_;
        label myProcNo_;
        label neighbProcNo_;


    // Private member functions

        //- Dissallow construction as copy
        procLduInterface(const procLduInterface&);


public:

    friend class LUscalarMatrix;


    // Constructors

        procLduInterface
        (
            const lduInterfaceField& interface,
            const scalarField& coeffs
        );

        procLduInterface(Istream& is);

        autoPtr<procLduInterface> clone()
        {
            notImplemented("procLduInterface::clone()");
            return autoPtr<procLduInterface>(NULL);
        }

        static autoPtr<procLduInterface> New(Istream& is)
        {
            return autoPtr<procLduInterface>(new procLduInterface(is));
        }


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const procLduInterface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
