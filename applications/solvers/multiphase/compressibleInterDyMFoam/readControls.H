#include "readTimeControls.H"

nAlphaSubCycles = readLabel(pimple.dict().lookup("nAlphaSubCycles"));

if (nAlphaSubCycles > 1 && pimple.nCorrPIMPLE() != 1)
{
    FatalErrorIn(args.executable())
        << "Sub-cycling alpha is only allowed for PISO, "
           "i.e. when the number of outer-correctors = 1"
        << exit(FatalError);
}

correctPhi = pimple.dict().lookupOrDefault("correctPhi", true);

checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);
