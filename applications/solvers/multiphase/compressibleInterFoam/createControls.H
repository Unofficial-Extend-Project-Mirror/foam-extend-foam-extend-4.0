#include "createTimeControls.H"

label nAlphaSubCycles
(
    readLabel(pimple.dict().lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1 && pimple.nCorrPIMPLE() != 1)
{
    FatalErrorIn(args.executable())
        << "Sub-cycling alpha is only allowed for PISO, "
           "i.e. when the number of outer-correctors = 1"
        << exit(FatalError);
}
