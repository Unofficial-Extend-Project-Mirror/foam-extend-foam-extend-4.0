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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "leastSquareInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquareInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        leastSquareInterpolation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::leastSquareInterpolation::calcAddressing() const
{
    if (addressingPtr_ || weightsPtr_)
    {
        FatalErrorIn
        (
            "void leastSquareInterpolation::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    const labelList& ac = overset().acceptorCells();

    addressingPtr_ = new labelListList(ac.size());
    labelListList& addr = *addressingPtr_;

    weightsPtr_ = new FieldField<Field, scalar>(ac.size());
    FieldField<Field, scalar>& w = *weightsPtr_;

    const labelList& donors = overset().localDonors();

    // Prepare a mask excluding hole and acceptor cells from the selection
    const labelList& hc = overset().holeCells();

    boolList donorMask(overset().mesh().nCells(), true);

    forAll (hc, hcI)
    {
        donorMask[hc[hcI]] = false;
    }

    forAll (ac, acI)
    {
        donorMask[ac[acI]] = false;
    }

    // Get cell-cell addressing
    const labelListList& cc = overset().mesh().cellCells();
    const vectorField& cellCentres = overset().mesh().cellCentres();

    forAll (addr, addrI)
    {
        // Start from donor cells
        DynamicList<label> curDonors(polyMesh::facesPerCell_);

        curDonors.append(donors[addrI]);

        // Consider neighbours
        const labelList& curNbrs = cc[donors[addrI]];

        forAll (curNbrs, nbrI)
        {
            // If neighbour is eligible, use it
            if (donorMask[curNbrs[nbrI]])
            {
                curDonors.append(curNbrs[nbrI]);
            }
        }
        Info<< "Acceptor: " << ac[addrI] << " root donor: " << donors[addrI]
            << " donors: " << curDonors << endl;
//         if (curDonors.size() > 1)
//         {
//             addr[addrI].transfer(curDonors.shrink());
//             const labelList& curAddr = addr[addrI];

//             // Calculate weights
//             w.set(addrI, new scalarField(curAddr.size(), 0));
//             scalarField& curW = w[addrI];

//             const vector& acCentre = cellCentres[ac[addrI]];

//             forAll (curAddr, caI)
//             {
//                 curW[caI] = 1/(mag(acCentre - cellCentres[curAddr[caI]]) + SMALL);
//             }

//             // Renormalise weights
//             curW /= gSum(curW);
//             Info<< "curW = " << sum(curW) << endl;
//         }
//         else
        {
            // Single donor addressing
            addr[addrI].setSize(1);
            addr[addrI][0] = donors[addrI];

            w.set(addrI, new scalarField(1, scalar(1)));
        }
    }
}


void Foam::leastSquareInterpolation::clearAddressing() const
{
    deleteDemandDrivenData(addressingPtr_);
    deleteDemandDrivenData(weightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::leastSquareInterpolation::leastSquareInterpolation
(
    const oversetMesh& overset,
    const dictionary& dict
)
:
    oversetInterpolation(overset, dict),
    addressingPtr_(NULL),
    weightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::leastSquareInterpolation::~leastSquareInterpolation()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::leastSquareInterpolation::localDonors() const
{
    return overset().localDonors();
}


const Foam::labelList& Foam::leastSquareInterpolation::remoteDonors() const
{
    return overset().remoteDonors();
}


const Foam::labelListList& Foam::leastSquareInterpolation::addressing() const
{
    if (!addressingPtr_)
    {
        calcAddressing();
    }

    return *addressingPtr_;
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::leastSquareInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


void Foam::leastSquareInterpolation::update()
{
    Info<< "leastSquareInterpolation::update()" << endl;
}


// ************************************************************************* //
