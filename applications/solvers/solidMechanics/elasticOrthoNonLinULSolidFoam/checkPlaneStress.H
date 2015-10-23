if(rheology.planeStress())
{
    //- add higher order terms
    volScalarField higherTerms = -0.5*volTensorField(gradDU&gradDU.T()).component(tensor::ZZ);

    forAll(gradDU.internalField(), celli)
    {
        gradDU.internalField()[celli][tensor::ZZ] =
            (
                (-C.internalField()[celli][symmTensor4thOrder::XXZZ]*DEpsilon.internalField()[celli][symmTensor::XX]
              - C.internalField()[celli][symmTensor4thOrder::YYZZ]*DEpsilon.internalField()[celli][symmTensor::YY]
            )
            /
            C.internalField()[celli][symmTensor4thOrder::ZZZZ])
          - higherTerms.internalField()[celli];
    }
    forAll(gradDU.boundaryField(), patchi)
    {
        forAll(gradDU.boundaryField()[patchi], facei)
        {
            gradDU.boundaryField()[patchi][facei][tensor::ZZ] =
                (
                    (
                        - C.boundaryField()[patchi][facei][symmTensor4thOrder::XXZZ]*
                          DEpsilon.boundaryField()[patchi][facei][symmTensor::XX]
                        - C.boundaryField()[patchi][facei][symmTensor4thOrder::YYZZ]*
                          DEpsilon.boundaryField()[patchi][facei][symmTensor::YY]
                    )
                    /
                    C.boundaryField()[patchi][facei][symmTensor4thOrder::ZZZZ]
                )
              - higherTerms.boundaryField()[patchi][facei];
        }
    }
}
