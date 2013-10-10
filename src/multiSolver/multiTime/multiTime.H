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
    Foam::multiTime

Description
    Gives read access to time directory format and precision.

SourceFiles
    multiTime.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef multiTime_H
#define multiTime_H

#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiTime Declaration
\*---------------------------------------------------------------------------*/

class multiTime
:
    public Time
{

public:

    TypeName("multiTime");


    // Constructors

        //- Construct given name, rootPath and casePath
        multiTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given dictionary, rootPath and casePath
        multiTime
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given endTime, rootPath and casePath
        multiTime
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


    // Destructor
    virtual ~multiTime();


    // Member functions
    
        // Access
        
            //- Time directory precision
            inline const int& precision() const
            {
                return Time::precision_;
            }
            
            //- Time directory format
            inline const Time::fmtflags& format() const
            {
                return Time::format_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
