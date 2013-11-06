word couplingScheme;

if (couplingProperties.found("couplingScheme"))
{
    couplingScheme = word(couplingProperties.lookup("couplingScheme"));

    if 
    (
        (couplingScheme == "IQN-ILS") 
     || (couplingScheme == "Aitken")
     || (couplingScheme == "FixedRelaxation")
    )
    {
        Info<< "Selecting coupling scheme " << couplingScheme << endl;
    }
    else
    {
        FatalErrorIn
        (
            "readFsiProperties"
        )   << "couplingScheme: " << couplingScheme
            << " is not a valid choice. "
            << "Options are: IQN-ILS, Aitken, FixedRelaxation"
            << abort(FatalError);
    }
}

// scalar couplingReuse = scalar
// (
//     readScalar
//     (
//         couplingProperties.lookup("couplingReuse")
//     )
// ); 

scalar interfaceDeformationLimit
(
    readScalar(couplingProperties.lookup("interfaceDeformationLimit"))
);


IOdictionary dynamicMeshDict
(
    IOobject
    (
        "dynamicMeshDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

// Switch twoDMotion(dynamicMeshDict.lookup("twoDMotion"));

scalar fsiRelaxationFactor;
scalar fsiRelaxationFactorMin;

fsiRelaxationFactorMin = scalar
(
    readScalar(couplingProperties.lookup("fsiRelaxationFactor"))
);
fsiRelaxationFactor = fsiRelaxationFactorMin;

scalar outerCorrTolerance
(
    readScalar(couplingProperties.lookup("outerCorrTolerance"))
); 

Switch fsi(couplingProperties.lookup("fsi"));
