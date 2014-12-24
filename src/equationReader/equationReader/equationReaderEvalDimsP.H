
    // Private member data

    // evaluateDims has an extra set of pointer functions that evaluateScalar
    // does not.  If equation::changeDimensions() is true, it is pointless to
    // perform a full evaluation, as the result will always be equation::
    // overrideDimensions().  Therefore evaluateDims has a pointer function
    // that either calls the standard evaluateDimsEnabled, or the shortcut
    // function evaluateDimsDisabled.  Since each equation may require
    // different treatment, there is a pointer function for each equation.  The
    // pointers are therefore in a list.

        //- Typedef the evaluateDims function pointer to make it available to
        //  PtrList
        typedef dimensionSet (Foam::equationReader::*evaluateDimsFunction)
        (
            const label equationIndex,
            const label maxStoreIndex
        ) const;

        // List of evaluateDims function pointers, indexed by equation
        mutable PtrList<evaluateDimsFunction> evaluateDimsFunctions_;

    // Private member functions

        //- Called when changeDimensions is false
        dimensionSet evaluateDimsEnabled
        (
            const label equationIndex,
            const label maxStoreIndex
        ) const;

        //- Called when changeDimensions is true
        dimensionSet evaluateDimsDisabled
        (
            const label equationIndex,
            const label maxStoreIndex
        ) const;


    // dimensionSet evaluation pointer functions

    void evalDimsNone
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsRetrieve
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsStore
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPlus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPlusDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMinus
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMinusDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsTimes
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsDivide
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPow
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPowDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSign
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSignDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsNeg
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMag
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLimit
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMinMod
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSqrtSumSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSqr
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPow3
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPow4
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPow5
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsPow6
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsInv
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSqrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsCbrt
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsHypot
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsHypotDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsExp
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsExpDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLog
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLogDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLog10
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLog10DimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSinDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsCos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsCosDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsTan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsTanDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAsin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAsinDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAcos
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAcosDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAtan
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAtanDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAtan2
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsSinhDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsCosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsCoshDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsTanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsTanhDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAsinh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAsinhDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAcosh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAcoshDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAtanh
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsAtanhDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsErf
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsErfDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsErfc
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsErfcDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLgamma
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsLgammaDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJ0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJ0DimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJ1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJ1DimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsJnDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsY0
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsY0DimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsY1
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsY1DimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsYn
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsYnDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMax
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMaxDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMin
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsMinDimCheck
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;

    void evalDimsStabilise
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storeIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;
