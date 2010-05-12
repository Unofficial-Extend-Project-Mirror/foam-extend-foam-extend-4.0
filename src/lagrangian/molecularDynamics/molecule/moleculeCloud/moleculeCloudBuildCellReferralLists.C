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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::moleculeCloud::buildCellReferralLists()
{
    Info<< nl << "Determining molecule referring schedule" << endl;

    const referredCellList& refIntL(referredInteractionList());

    DynamicList<label> referralProcs;

    // Run through all referredCells to build list of interacting processors

    forAll(refIntL, rIL)
    {
        const referredCell& rC(refIntL[rIL]);

        if (findIndex(referralProcs, rC.sourceProc()) == -1)
        {
            referralProcs.append(rC.sourceProc());
        }
    }

    referralProcs.shrink();

//     Pout << "referralProcs: " << nl << referralProcs << endl;

    List<DynamicList<label> > cellSendingReferralLists(referralProcs.size());

    List<DynamicList<DynamicList<label> > >
        cellReceivingReferralLists(referralProcs.size());

    // Run through all referredCells again building up send and receive info

    forAll(refIntL, rIL)
    {
        const referredCell& rC(refIntL[rIL]);

        label rPI = findIndex(referralProcs, rC.sourceProc());

        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        label existingSource = findIndex(sRL, rC.sourceCell());

        // Check to see if this source cell has already been allocated to
        // come to this processor.  If not, add the source cell to the sending
        // list and add the current referred cell to the receiving list.

        // It shouldn't be possible for the sending and receiving lists to be
        // different lengths, because their append operations happen at the
        // same time.

        if (existingSource == -1)
        {
            sRL.append(rC.sourceCell());

            rRL.append
            (
                DynamicList<label> (labelList(1,rIL))
            );
        }
        else
        {
            rRL[existingSource].append(rIL);

            rRL[existingSource].shrink();
        }
    }

    forAll(referralProcs, rPI)
    {
        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        sRL.shrink();

        rRL.shrink();
    }

    // It is assumed that cell exchange is reciprocal, if proc A has cells to
    // send to proc B, then proc B must have some to send to proc A.

    cellReceivingReferralLists_.setSize(referralProcs.size());

    cellSendingReferralLists_.setSize(referralProcs.size());

    forAll(referralProcs, rPI)
    {
        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        labelListList translLL(rRL.size());

        forAll(rRL, rRLI)
        {
            translLL[rRLI] = rRL[rRLI];
        }

        cellReceivingReferralLists_[rPI] = receivingReferralList
        (
            referralProcs[rPI],
            translLL
        );
    }

    // Send sendingReferralLists to each interacting processor.

    forAll(referralProcs, rPI)
    {

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                OPstream toInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                toInteractingProc << sendingReferralList
                (
                    Pstream::myProcNo(),
                    sRL
                );
            }
        }
    }

    // Receive sendingReferralLists from each interacting processor.

    forAll(referralProcs, rPI)
    {
        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                IPstream fromInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                fromInteractingProc >> cellSendingReferralLists_[rPI];
            }
        }
        else
        {
            cellSendingReferralLists_[rPI] = sendingReferralList
            (
                Pstream::myProcNo(),
                cellSendingReferralLists[rPI]
            );
        }
    }

//     Pout << "receiving list: " << nl << cellReceivingReferralLists_ << endl;

//     Pout << "sending list: " << nl << cellSendingReferralLists_ << endl;
}


// ************************************************************************* //
