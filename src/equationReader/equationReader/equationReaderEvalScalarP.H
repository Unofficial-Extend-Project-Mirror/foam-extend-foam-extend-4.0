    // scalar evaluation pointer functions

    void evalScalarNone
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarRetrieve
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarStore
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPlus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarMinus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarTimes
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarDivide
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPow
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSign
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarNeg
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarMag
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarLimit
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarMinMod
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSqrtSumSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPow3
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPow4
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPow5
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarPow6
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarInv
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSqrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarCbrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarHypot
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarExp
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarLog
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarLog10
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarCos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarTan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAsin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAcos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAtan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAtan2
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarSinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarCosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarTanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAsinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAcosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarAtanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarErf
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarErfc
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarLgamma
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarJ0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarJ1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarJn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarY0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarY1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarYn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarMax
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarMin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

    void evalScalarStabilise
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        scalar& x,
        scalar source
    ) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
