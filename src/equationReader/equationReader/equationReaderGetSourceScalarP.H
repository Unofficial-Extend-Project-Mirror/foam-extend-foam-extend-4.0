    // Get source scalar functions
    scalar getScalarSrcNone
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcStorage
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcActiveSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcEquation
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcEquationCircRefDetect
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcInternalScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcDictSourceDScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcDictSourceScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcScalarSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcScalarFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcVectorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcVectorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcDiagTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcDiagTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcSymmTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcSymmTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcSphericalTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    scalar getScalarSrcSphericalTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

