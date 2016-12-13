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

\*---------------------------------------------------------------------------*/

#include "procLduInterface.H"
#include "lduInterfaceField.H"
#include "cyclicLduInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procLduInterface::procLduInterface
(
    const lduInterfaceField& interface,
    const scalarField& coeffs
)
:
    faceCells_(interface.coupledInterface().faceCells()),
    coeffs_(coeffs),
    myProcNo_(-1),
    neighbProcNo_(-1)
{
    if (isA<processorLduInterface>(interface.coupledInterface()))
    {
        const processorLduInterface& pldui =
            refCast<const processorLduInterface>(interface.coupledInterface());

        myProcNo_ = pldui.myProcNo();
        neighbProcNo_ = pldui.neighbProcNo();
    }
    else if (isA<cyclicLduInterface>(interface.coupledInterface()))
    {
    }
    else
    {
        FatalErrorIn
        (
            "procLduInterface::procLduInterface"
            "(const lduInterfaceField&, const scalarField&"
        )   << "unknown lduInterface type "
                << interface.coupledInterface().type()
            << exit(FatalError);
    }
}


Foam::procLduInterface::procLduInterface(Istream& is)
:
    faceCells_(is),
    coeffs_(is),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is))
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const procLduInterface& cldui)
{
    os  << cldui.faceCells_
        << cldui.coeffs_
        << cldui.myProcNo_
        << cldui.neighbProcNo_;

    return os;
}


// ************************************************************************* //
