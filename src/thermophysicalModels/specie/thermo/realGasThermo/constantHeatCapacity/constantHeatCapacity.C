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

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#include "constantHeatCapacity.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::constantHeatCapacity<equationOfState>::constantHeatCapacity(Istream& is)
:
    equationOfState(is),
    //CL: reads specific perfect gas heat capacity
    Cp0_(readScalar(is)),
    cp0_(Cp0_*this->W()),
    //CL: values for some need terms at std
    e0_std(e0(this->Tstd())),
    s0_std(s0(this->Tstd())),
    integral_p_dv_std(this->integral_p_dv(this->rhostd(),this->Tstd())),
    integral_dpdT_dv_std(this->integral_dpdT_dv(this->rhostd(),this->Tstd())),
    //CL: cp @ STD (needed to limit cp for stability)
    cp_std(this->cp_nonLimited(this->rhostd(), this->Tstd()))
{
    is.check("constantHeatCapacity::constantHeatCapacity(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const constantHeatCapacity<equationOfState>& np
)
{
    os  << static_cast<const equationOfState&>(np) << tab
        << np.Cp0_;
    os.check("Ostream& operator<<(Ostream& os, const constantHeatCapacity& np)");
    return os;
}


// ************************************************************************* //
