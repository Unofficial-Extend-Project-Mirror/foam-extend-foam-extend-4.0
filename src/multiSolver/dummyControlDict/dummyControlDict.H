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
    Foam::dummyControlDict

Description
    This class is a workaround that allows a class to have its own independent
    objectRegistry.  That requires the creation of a new Time object, which
    requires a controlDict for construction.  dummyControlDict creates a dummy
    controlDict to feed the new Time object.  The only dictionary values it
    handles are deltaT, writeFrequency, and optionally, the Time class statics,
    format and precision.

SourceFiles
    dummyControlDict.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef dummyControlDict_H
#define dummyControlDict_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dummyControlDict Declaration
\*---------------------------------------------------------------------------*/

class dummyControlDict
:
    public dictionary
{

public:

    TypeName("dummyControlDict");


    // Constructors

        //- Construct, creating the bare minimum controlDict
        dummyControlDict();

        //- Construct, reading the multiControlDict data from a file
        explicit dummyControlDict(const fileName&);
        
        //- Construct given the multiControlDict dictionary
        explicit dummyControlDict(const dictionary&);


    // Destructor
    virtual ~dummyControlDict();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
