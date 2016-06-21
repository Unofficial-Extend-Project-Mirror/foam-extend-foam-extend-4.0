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
    coupledLduMatrix

Description
    Collection of lduMatrices solved together as a block system

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupledLduMatrix.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledLduMatrix, 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given size
Foam::coupledLduMatrix::coupledLduMatrix(const label size)
:
    PtrList<lduMatrix>(size)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledLduMatrix::~coupledLduMatrix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coupledLduMatrix::diagonal() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool diag = true;

    forAll (matrices, matrixI)
    {
        diag = diag && matrices[matrixI].diagonal();
    }

    return diag;
}


bool Foam::coupledLduMatrix::symmetric() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool sym = true;

    forAll (matrices, matrixI)
    {
        sym =
            (sym && matrices[matrixI].diagonal())
         || (sym && matrices[matrixI].symmetric());
    }

    return sym;
}


bool Foam::coupledLduMatrix::asymmetric() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool asym = false;

    forAll (matrices, matrixI)
    {
        asym = (asym || matrices[matrixI].asymmetric());
    }

    return asym;
}


void Foam::coupledLduMatrix::Amul
(
    FieldField<Field, scalar>& result,
    const FieldField<Field, scalar>& x,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    // Reset product to zero
    result = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        bouCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );

    forAll (matrices, rowI)
    {
        matrices[rowI].AmulCore(result[rowI], x[rowI]);
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        bouCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );
}


void Foam::coupledLduMatrix::Tmul
(
    FieldField<Field, scalar>& result,
    const FieldField<Field, scalar>& x,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    // Reset product to zero
    result = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        intCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );

    forAll (matrices, rowI)
    {
        matrices[rowI].TmulCore(result[rowI], x[rowI]);
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        intCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );
}


void Foam::coupledLduMatrix::initMatrixInterfaces
(
    const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const FieldField<Field, scalar>& x,
    FieldField<Field, scalar>& result,
    const direction cmpt,
    const bool switchToLhs
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    // Note.  The comms design requires all non-processor interfaces
    // to be updated first, followed by the update of processor
    // interfaces.  The reason is that non-processor coupled
    // interfaces require a complex comms pattern involving more than
    // pairwise communications.
    // Under normal circumstances this is achieved naturally, since
    // processor interfaces come last on the list and other coupled
    // interfaces execute complex comms at init() level.
    // For coupled matrices, the update loop needs to be split over
    // all matrices by hand
    // Bug fix: Zeljko Tukovic, 7/Apr/2015

    // Init update all non-processor coupled interfaces
    forAll (matrices, rowI)
    {
        if
        (
            Pstream::defaultCommsType() == Pstream::blocking
         || Pstream::defaultCommsType() == Pstream::nonBlocking
        )
        {
            forAll (interfaces[rowI], interfaceI)
            {
                if (interfaces[rowI].set(interfaceI))
                {
                    if
                    (
                       !isA<processorLduInterfaceField>
                        (
                            interfaces[rowI][interfaceI]
                        )
                    )
                    {
                        interfaces[rowI][interfaceI].initInterfaceMatrixUpdate
                        (
                            x[rowI],
                            result[rowI],
                            matrices[rowI],
                            coupleCoeffs[rowI][interfaceI],
                            cmpt,
                            static_cast<const Pstream::commsTypes>
                            (
                                Pstream::defaultCommsType()
                            ),
                            switchToLhs
                        );
                    }
                }
            }
        }
        else if (Pstream::defaultCommsType() == Pstream::scheduled)
        {
            // ERROR: Does not work with scheduled comms.
            // To investigate.  HJ, 11/Jun/2015
            forAll (matrices, rowI)
            {
                const lduSchedule& patchSchedule =
                    matrices[rowI].patchSchedule();

                // Loop over the "global" patches are on the list of
                // interfaces but beyond the end of the schedule
                // which only handles "normal" patches
                for
                (
                    label interfaceI = patchSchedule.size()/2;
                    interfaceI < interfaces.size();
                    interfaceI++
                )
                {
                    if (interfaces[rowI].set(interfaceI))
                    {
                        if
                        (
                            !isA<processorLduInterfaceField>
                            (
                                interfaces[rowI][interfaceI]
                            )
                        )
                        {
                            interfaces[rowI][interfaceI].
                                initInterfaceMatrixUpdate
                                (
                                    x[rowI],
                                    result[rowI],
                                    matrices[rowI],
                                    coupleCoeffs[rowI][interfaceI],
                                    cmpt,
                                    static_cast<const Pstream::commsTypes>
                                    (
                                        Pstream::defaultCommsType()
                                    ),
                                    switchToLhs
                                );
                        }
                    }
                }
            }
        }
        else
        {
            FatalErrorIn("void coupledLduMatrix::initMatrixInterfaces")
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType()]
                << exit(FatalError);
        }
    }

    // Init update for all processor interfaces
    forAll (matrices, rowI)
    {
        if
        (
            Pstream::defaultCommsType() == Pstream::blocking
         || Pstream::defaultCommsType() == Pstream::nonBlocking
        )
        {
            forAll (interfaces[rowI], interfaceI)
            {
                if (interfaces[rowI].set(interfaceI))
                {
                    if
                    (
                        isA<processorLduInterfaceField>
                        (
                            interfaces[rowI][interfaceI]
                        )
                    )
                    {
                        interfaces[rowI][interfaceI].initInterfaceMatrixUpdate
                        (
                            x[rowI],
                            result[rowI],
                            matrices[rowI],
                            coupleCoeffs[rowI][interfaceI],
                            cmpt,
                            static_cast<const Pstream::commsTypes>
                            (
                                Pstream::defaultCommsType()
                            ),
                            switchToLhs
                        );
                    }
                }
            }
        }
        else if (Pstream::defaultCommsType() == Pstream::scheduled)
        {
            // ERROR: Does not work with scheduled comms.
            // To investigate.  HJ, 11/Jun/2015
            forAll (matrices, rowI)
            {
                const lduSchedule& patchSchedule =
                    matrices[rowI].patchSchedule();

                // Loop over the "global" patches are on the list of
                // interfaces but beyond the end of the schedule
                // which only handles "normal" patches
                for
                (
                    label interfaceI = patchSchedule.size()/2;
                    interfaceI < interfaces.size();
                    interfaceI++
                )
                {
                    if (interfaces[rowI].set(interfaceI))
                    {
                        if
                        (
                            isA<processorLduInterfaceField>
                            (
                                interfaces[rowI][interfaceI]
                            )
                        )
                        {
                            interfaces[rowI][interfaceI].
                                initInterfaceMatrixUpdate
                                (
                                    x[rowI],
                                    result[rowI],
                                    matrices[rowI],
                                    coupleCoeffs[rowI][interfaceI],
                                    cmpt,
                                    static_cast<const Pstream::commsTypes>
                                    (
                                        Pstream::defaultCommsType()
                                    ),
                                    switchToLhs
                                );
                        }
                    }
                }
            }
        }
        else
        {
            FatalErrorIn("void coupledLduMatrix::initMatrixInterfaces")
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType()]
                << exit(FatalError);
        }
    }
}


void Foam::coupledLduMatrix::updateMatrixInterfaces
(
    const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const FieldField<Field, scalar>& x,
    FieldField<Field, scalar>& result,
    const direction cmpt,
    const bool switchToLhs
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    forAll (matrices, rowI)
    {
        matrices[rowI].updateMatrixInterfaces
        (
            coupleCoeffs[rowI],
            interfaces[rowI],
            x[rowI],
            result[rowI],
            cmpt,
            switchToLhs
        );
    }
}


// ************************************************************************* //
