scalar xMax = 0;

scalar yMax = 0;

scalar zMax = 0;

scalar xMin = 0;

scalar yMin = 0;

scalar zMin = 0;

label xMaxPtLabel = 0;

label yMaxPtLabel = 0;

label zMaxPtLabel = 0;

label xMinPtLabel = 0;

label yMinPtLabel = 0;

label zMinPtLabel = 0;

forAll (cellZoneI[cZ], nC)
{
    const labelList& cellPointsJ = mesh_.cellPoints()[cellZoneI[cZ][nC]];

    forAll(cellPointsJ, nP)
    {
        const point& ptI = mesh_.points()[cellPointsJ[nP]];

        const label& ptILabel = cellPointsJ[nP];

        if (ptI.x() > xMax || nC == 0)
        {
            xMax = ptI.x();
            xMaxPtLabel = ptILabel;
        }
        if (ptI.y() > yMax || nC == 0)
        {
            yMax = ptI.y();
            yMaxPtLabel = ptILabel;
        }
        if (ptI.z() > zMax || nC == 0)
        {
            zMax = ptI.z();
            zMaxPtLabel = ptILabel;
        }
        if (ptI.x() < xMin || nC == 0)
        {
            xMin = ptI.x();
            xMinPtLabel = ptILabel;
        }
        if (ptI.y() < yMin || nC == 0)
        {
            yMin = ptI.y();
            yMinPtLabel = ptILabel;
        }
        if (ptI.z() < zMin || nC == 0)
        {
            zMin = ptI.z();
            zMinPtLabel = ptILabel;
        }
    }
}

// Info << "Xmax: label = " << xMaxPtLabel2 << "; vector = " <<mesh_.points()[xMaxPtLabel2]
//      <<"; x-component = " << mesh_.points()[xMaxPtLabel2].x() << endl;
// Info << "Ymax: label = " << yMaxPtLabel2 << "; vector = " <<mesh_.points()[yMaxPtLabel2]
//      <<"; y-component = " << mesh_.points()[yMaxPtLabel2].y() << endl;
// Info << "Zmax: label = " << zMaxPtLabel2 << "; vector = " <<mesh_.points()[zMaxPtLabel2]
//      <<"; z-component = " << mesh_.points()[zMaxPtLabel2].z() << endl;
//
// Info << "Xmin: label = " << xMinPtLabel << "; vector = " <<mesh_.points()[xMinPtLabel]
//      <<"; x-component = " << mesh_.points()[xMinPtLabel].x() << endl;
// Info << "Ymin: label = " << yMinPtLabel << "; vector = " <<mesh_.points()[yMinPtLabel]
//      <<"; y-component = " << mesh_.points()[yMinPtLabel].y() << endl;
// Info << "Zmin: label = " << zMinPtLabel << "; vector = " <<mesh_.points()[zMinPtLabel]
//      <<"; z-component = " << mesh_.points()[zMinPtLabel].z() << endl;

scalar xMid =
    (mesh_.points()[xMaxPtLabel].x()
  + mesh_.points()[xMinPtLabel].x()) / 2;

scalar yMid =
    (mesh_.points()[yMaxPtLabel].y()
  + mesh_.points()[yMinPtLabel].y()) / 2;

scalar zMid =
    (mesh_.points()[zMaxPtLabel].z()
  + mesh_.points()[zMinPtLabel].z()) / 2;

vector rS(xMid, yMid, zMid);

// Info << "\t The Estimated Starting Point: " << rS << endl;

