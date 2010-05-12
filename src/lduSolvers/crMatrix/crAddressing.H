/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    crAddressing

Description
    Compressed row matrix addressing

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    crAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef crAddressing_H
#define crAddressing_H

#include "labelList.H"
#include "refCount.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class crAddressing Declaration
\*---------------------------------------------------------------------------*/

class crAddressing
:
    public refCount
{
    // Private data

        //- Number of rows
        label nRows_;

        //- Number of columns
        label nCols_;

        //- Row array.
        //  Provides start index for each row, dimensioned to nRows + 1
        labelList row_;

        //- Column array
        labelList col_;

    // Private Member Functions

        //- Set row count
        void setRowCount(const labelList& count);


public:

    // Constructors

        //- Construct given size.  Column and coefficients set later
        crAddressing
        (
            const label nRows,
            const label nCols,
            const labelList& count
        );

        //- Construct from components
        crAddressing
        (
            const label nRows,
            const label nCols,
            const labelList& row,
            const labelList& col
        );

        //- Construct as copy
        crAddressing(const crAddressing&);

        //- Construct from Istream
        crAddressing(Istream&);


    // Destructor - default


    // Member Functions

        // Access

            //- Return number of rows
            label nRows() const
            {
                return nRows_;
            }

            //- Return number of columns
            label nCols() const
            {
                return nCols_;
            }

            //- Return number of coefficients
            label nEntries() const
            {
                return col_.size();
            }

            //- Return row array
            const labelList& row() const
            {
                return row_;
            }

            //- Return column array
            const labelList& col() const
            {
                return col_;
            }


        // Edit

            //- Return column array
            labelList& col()
            {
                return col_;
            }


        // Operations

            //- Return transpose addressing
            tmp<crAddressing> T() const;


    // Member Operators

        void operator=(const crAddressing&);

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const crAddressing&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
