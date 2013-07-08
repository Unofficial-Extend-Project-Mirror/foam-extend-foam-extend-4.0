//****************************************************//
// Read the contact tol and max contact correctors
//****************************************************//
const IOdictionary& contactControl =
  IOobject
  (
   "contactProperties",
   U.time().constant(),
   U.db(),
   IOobject::MUST_READ,
   IOobject::NO_WRITE
   );


//- The frequnecy the conatct is corrected inside the momentum loop
int uEqnContactCorrFreq(readInt(contactControl.lookup("innerContactCorrectFreq")));
