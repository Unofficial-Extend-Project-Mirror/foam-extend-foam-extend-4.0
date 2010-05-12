bool meshMoving = true;

if (Times.size() > 2)
{
   for(label n2=2; n2<Times.size(); n2++)
   {
       IOobject tmpPoints
       (
           "points",
           Times[n2].name(),
           polyMesh::meshSubDir,
           mesh,
           IOobject::NO_READ
       );
       meshMoving = meshMoving && tmpPoints.headerOk();
   }
}
else
{
    meshMoving = false;
}
