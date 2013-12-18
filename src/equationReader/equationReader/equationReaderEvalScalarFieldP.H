    // scalar evaluation pointer functions

    void evalScalarFieldNone
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldRetrieve
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldStore
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPlus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldMinus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldTimes
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldDivide
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPow
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSign
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldNeg
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldMag
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldLimit
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldMinMod
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSqrtSumSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPow3
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPow4
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPow5
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldPow6
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldInv
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSqrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldCbrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldHypot
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldExp
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldLog
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldLog10
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldCos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldTan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAsin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAcos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAtan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAtan2
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldSinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldCosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldTanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAsinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAcosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldAtanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldErf
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldErfc
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldLgamma
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldJ0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldJ1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldJn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldY0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldY1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldYn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldMax
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldMin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    void evalScalarFieldStabilise
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalarField& x,
        const scalarField& source
    ) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
