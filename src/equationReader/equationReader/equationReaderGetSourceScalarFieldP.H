    // Get source scalar functions
    const scalarField& getScalarFieldSrcNone
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcStorage
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcActiveSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcEquation
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcEquationCircRefDetect
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcInternalScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcDictSourceDScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcDictSourceScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcScalarSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcScalarFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcVectorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcVectorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcDiagTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcDiagTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcSymmTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcSymmTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcSphericalTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    const scalarField& getScalarFieldSrcSphericalTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

