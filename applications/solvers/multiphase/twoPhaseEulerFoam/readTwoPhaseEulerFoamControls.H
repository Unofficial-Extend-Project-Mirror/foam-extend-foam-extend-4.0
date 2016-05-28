#   include "readTimeControls.H"

    int nAlphaCorr(readInt(piso.dict().lookup("nAlphaCorr")));

    Switch correctAlpha(piso.dict().lookup("correctAlpha"));
