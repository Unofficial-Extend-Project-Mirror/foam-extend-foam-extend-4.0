    // Motion is a vectorField of all moving boundary points
    vectorField motion(ms.movingPoints().size(), vector::zero);

    vectorField oldPoints = ms.movingPoints();

    scalar oldTime = time_.value() - time_.deltaT().value();
    scalar curTime = time_.value();

    scalar alphaOld   = 0.0;
    scalar alphaCur   = 0.0;

//     scalar rotationAmplitude_(0.707);
//     scalar rotationFrequency_(0.0);
//     vector translationFrequency_(0.0,0.25,0.0);
//     vector translationAmplitude_(0.0,5.0,0.0);
//     vector initialRotationOrigin_(0.0,0.0,0.0);

    scalar pi=3.141592;

    alphaOld = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*oldTime);
    alphaCur = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*curTime);

    vector translationVector
    (
        translationAmplitude_[0]*
        (
            Foam::sin(2*pi*translationFrequency_[0]*curTime)
          - Foam::sin(2*pi*translationFrequency_[0]*oldTime)
        ),
        translationAmplitude_[1]*
        (
            Foam::sin(2*pi*translationFrequency_[1]*curTime)
          - Foam::sin(2*pi*translationFrequency_[1]*oldTime)
        ),
        0
    );

    tensor RzOld
    (
        Foam::cos(alphaOld), -Foam::sin(alphaOld), 0,
        Foam::sin(alphaOld),  Foam::cos(alphaOld), 0,
        0, 0, 1
    );

    tensor RzCur
    (
        Foam::cos(alphaCur), -Foam::sin(alphaCur), 0,
        Foam::sin(alphaCur),  Foam::cos(alphaCur), 0,
        0, 0, 1
    );

    vectorField rotationField
    (
        (RzCur - RzOld)
      & (statPoints_ - initialRotationOrigin_)
    );

    motion = translationVector + rotationField;
