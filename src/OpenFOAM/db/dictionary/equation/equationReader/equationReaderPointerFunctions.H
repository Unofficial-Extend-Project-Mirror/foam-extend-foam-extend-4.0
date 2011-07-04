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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Get source pointer functions
    dimensionedScalar getSourceNone
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceDictSourceDScalar
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceDictSourceScalar
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceExternalDScalar
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceExternalScalar
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceExternalScalarList
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceInternalScalar
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceEquation
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceEquationCircRefDetect
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );

    dimensionedScalar getSourceStorage
    (
        equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    );
    

    // Evaluation pointer functions
    void evalNone
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalRetrieve
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalRetrieveChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalStore
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPlus
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPlusDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPlusChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinus
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinusDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinusChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTimes
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalDivide
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPow
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPowDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPowChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSign
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPos
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalNeg
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMag
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLimit
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinMod
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSqrtSumSqr
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSqr
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPow3
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPow4
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalPow5
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalInv
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSqrt
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCbrt
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalHypot
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalHypotDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalHypotChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalExp
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalExpDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalExpChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLog
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLogDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLogChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLog10
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLog10DimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLog10ChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSin
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSinDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSinChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCos
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCosDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCosChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTan
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTanDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTanChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsin
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsinDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsinChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcos
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcosDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcosChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtan
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtanDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtanChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSinh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSinhDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalSinhChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCosh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCoshDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalCoshChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTanh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTanhDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalTanhChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsinh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsinhDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAsinhChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcosh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcoshDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAcoshChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtanh
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtanhDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalAtanhChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErf
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErfDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErfChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErfc
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErfcDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalErfcChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLgamma
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLgammaDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalLgammaChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ0
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ0DimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ0ChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ1
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ1DimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJ1ChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJn
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJnDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalJnChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY0
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY0DimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY0ChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY1
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY1DimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalY1ChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalYn
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalYnDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalYnChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMax
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMaxDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMaxChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMin
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinDimCheck
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalMinChangeDimensions
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

    void evalStabilise
    (
        equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionedScalar& ds,
        dimensionedScalar source
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
