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
    Foam::ash

Description
    Coal ash

SourceFiles
    ash.C

\*---------------------------------------------------------------------------*/

#ifndef solid_ash_H
#define solid_ash_H

#include "solid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class ash Declaration
\*---------------------------------------------------------------------------*/

class ash
:
    public solid
{

public:

    //- Runtime type information
    TypeName("ash");


    // Constructors

        //- Construct null
        ash()
        :
//   rho       Cp          K        Hf      e
// [kg/m3] [J/(kg.K)] [W/(m.K)]   [J/kg]   [-]
            solid(2010, 710, 0.04, 0.0, 1.0) // TODO
        {
            WarningIn("ash::ash()")
                << "Properties of ash need to be checked!!!"
                << endl;
        }

        ash
        (
            const solid& s
        )
        :
            solid(s)
        {}

        //- Construct from Istream
        ash(Istream& is)
        :
            solid(is)
        {}

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            solid::writeData(os);
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const ash& s)
        {
            s.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
