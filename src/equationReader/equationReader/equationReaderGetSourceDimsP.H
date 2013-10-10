    // Get source dimension functions
    dimensionSet getDimsSrcNone
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcStorage
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcActiveSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcEquation
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcEquationCircRefDetect
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcInternalScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcDictSourceDScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcDictSourceScalar
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcScalarSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcScalarFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcVectorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcVectorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcDiagTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcDiagTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcSymmTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcSymmTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcSphericalTensorSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    dimensionSet getDimsSrcSphericalTensorFieldSource
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

