Switch solidInterfaceCorr(false);

solidInterface* solidInterfacePtr(NULL);
IOReferencer<solidInterface>* solidInterfaceRefPtr(NULL);

{
    const dictionary& stressControl =
      mesh.solutionDict().subDict("stressedFoam");

    solidInterfaceCorr = Switch(stressControl.lookup("solidInterface"));

    if(solidInterfaceCorr)
      {
	Info << "Creating solid interface correction" << endl;
	solidInterfacePtr = new solidInterface(mesh, rheology);
	solidInterfacePtr->modifyProperties(muf, lambdaf);

	// leastSquaresSolidInterfaceGrad needs to lookup solidInterface object
	// from the object registry so I will use IOReferencer to register
	// the solidInterface object
	solidInterface& solidInterfaceObj = *solidInterfacePtr;
	solidInterfaceRefPtr =
	  new IOReferencer<solidInterface>
	  (
	   IOobject
	   (
	    "solidInterface",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,  /*must be NO_READ*/
	    IOobject::NO_WRITE  /*must be NO_WRITE*/
	    ),
	   solidInterfaceObj
	   );

	//- solidInterface needs muf and lambdaf to be used for divSigmaExp
	if(divSigmaExpMethod != "surface" && divSigmaExpMethod != "decompose")
	  {
	    FatalError << "divSigmaExp must be decompose or surface when solidInterface is on"
		       << exit(FatalError);
	  }

	// check grad scheme
	if(word(mesh.gradSchemes().lookup("grad(U)")) != "leastSquaresSolidInterface")
	  {
	    Warning << "The grad(U) gradScheme should be leastSquaresSolidInterface "
		    << "for the solidInterface procedure to be correct"
		    << endl;
	  }
      }
}

//IOReferencer<solidInterface>& solidInterfaceRef = *solidInterfaceRefPtr;;
