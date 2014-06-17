// dirichletNeumannCalculateSlaveFaceNormals.H

   if (nonLinear == "off")
     {
       // undeformed normals
       // remember that the mesh has not moved,
         // only the global face zone patches
       // are moved to the deformed position

       // slaveFaceNormals =
         //-1*masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
       //        (
       //         mesh.boundaryMesh()[masterPatchIndex].faceNormals()
       //         );
       // slaveFaceNormals /= mag(slaveFaceNormals);

       // OPTION 2
       // Use deformed slave face normals
       // Get them from global face zone patches
       // const vectorField& globalSlaveFaceNormals =
         // slaveFaceZonePatch.faceNormals();
       // const label slavePatchStart
       //    = mesh.boundaryMesh()[slavePatchIndex].start();
       // forAll(slaveFaceNormals, facei)
       //    {
       //      slaveFaceNormals[facei] =
       //        globalSlaveFaceNormals
       //        [
       //         mesh.faceZones()[slaveFaceZoneID()].whichFace(
         //slavePatchStart + facei
         //  )
       //         ];
       //    }
       // // make sure they are unity
       // slaveFaceNormals /= mag(slaveFaceNormals);

       // OPTION 3
       // globally interpolate master normals to the slave
       // then get the local normals

       // undeformed master normals
       // const vectorField& actualSlaveFaceNormals =
         // mesh.boundaryMesh()[slavePatchIndex].faceNormals();
       // vectorField globalMasterFaceNormals
         //(masterFaceZonePatch.size(), vector::zero);
       // const label masterPatchStart
       //    = mesh.boundaryMesh()[masterPatchIndex].start();
       // // master undeformed face normals
       // vectorField masterFaceNormals =
         // mesh.boundaryMesh()[masterPatchIndex].faceNormals();
       // // put field into global
       // forAll(masterFaceNormals, i)
       //    {
       //      globalMasterFaceNormals
         // [
         // mesh.faceZones()[masterFaceZoneID()
         // ].whichFace(masterPatchStart + i)] =
       //        masterFaceNormals[i];
       //    }

       // //- exchange parallel data
       // // sum because each face is only on one proc
       // reduce(globalMasterFaceNormals, sumOp<vectorField>());

       // use deformed master normals
       const vectorField& globalMasterFaceNormals =
           masterFaceZonePatch.faceNormals();

       // interpolate master normals to the slave and reverse the direction
       // we can use inverseDistance or GGI for the interpolation
       // GGI is much better is we want the actual master normals
       vectorField globalSlaveFaceNormals
           (slaveFaceZonePatch.size(), vector::zero);
       if (ggiInterpolatorPtr)
     {
       //Info << "interpolating master normals with GGI" << endl;
       globalSlaveFaceNormals =
         -ggiInterpolatorPtr->masterToSlave(globalMasterFaceNormals);
     }
       else // inverse distance
     {
       globalSlaveFaceNormals =
         -masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
         (
          globalMasterFaceNormals
          );
     }
       vectorField actualSlaveFaceNormals =
     slaveFaceZonePatch.faceNormals();

       // put global back into local
       const label slavePatchStart
     = mesh.boundaryMesh()[slavePatchIndex].start();
       forAll(slaveFaceNormals, facei)
             {
               slaveFaceNormals[facei] =
                 globalSlaveFaceNormals
                 [
                  mesh.faceZones()[slaveFaceZoneID()].whichFace
                  (slavePatchStart + facei)
                  ];

       if (mag(slaveFaceNormals[facei]) < SMALL)
         {
           // interpolation sometimes does not work for far
           // away faces so we will use actual slave normals
           // but doesn't really matter as the face is probably
           // not in contact
           slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
         }
       else
         {
           // make sure they are unity
           slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
         }
             }
     }
     else if
     (
         nonLinear == "updatedLagrangian"
      || nonLinear == "totalLagrangian"
     )
     {
       // a few different ways to calculate deformed normals

       // OPTION 1
       // use deformed slave normals calculated using F
       // calculate deformed normals using deformation gradient
       // const volTensorField& gradField =
       //    mesh.objectRegistry::lookupObject<volTensorField>
         //("grad("+fieldName+")");
       //        // deformation gradient
       //        tensorField F = gradField.boundaryField()[slavePatchIndex];
       //        if (fieldName == "DU")
       //    {
       //      const volTensorField& gradU =
       //        mesh.objectRegistry::lookupObject<volTensorField>("grad(U)");
       //      F += gradU.boundaryField()[slavePatchIndex];
       //    }
       //        const tensorField Finv = inv(I + F);
       //        // deformed normals
       //        slaveFaceNormals =
         //Finv & mesh.boundaryMesh()[slavePatchIndex].faceNormals();
       //        slaveFaceNormals /= mag(slaveFaceNormals);

       // OPTION 2
       // Use deformed slave face normals
       // Get them from global face zone patches
       // const vectorField& globalSlaveFaceNormals =
         //slaveFaceZonePatch.faceNormals();
       // const label slavePatchStart
       //    = mesh.boundaryMesh()[slavePatchIndex].start();
       // forAll(slaveFaceNormals, facei)
       //    {
       //      slaveFaceNormals[facei] =
       //        globalSlaveFaceNormals
       //        [
       //         mesh.faceZones()[slaveFaceZoneID()].whichFace
         //(slavePatchStart + facei)
       //         ];
       //    }
       // // make sure they are unity
       // slaveFaceNormals /= mag(slaveFaceNormals);

       // OPTION 3
       // Use deformed master normals
       // globally interpolate deformed master normals to the slave
       // then get the local normals
       // globally interpolate master normals to the slave
       // then get the local normals

       // undeformed master normals
       //const vectorField& actualSlaveFaceNormals =
         // mesh.boundaryMesh()[slavePatchIndex].faceNormals();
       // vectorField globalMasterFaceNormals
         // (masterFaceZonePatch.size(), vector::zero);
       // const label masterPatchStart
       //    = mesh.boundaryMesh()[masterPatchIndex].start();
       // // master undeformed face normals
       // vectorField masterFaceNormals =
         // mesh.boundaryMesh()[masterPatchIndex].faceNormals();
       // // put field into global
       // forAll(masterFaceNormals, i)
       //    {
       //      globalMasterFaceNormals
         // [
         // mesh.faceZones()[masterFaceZoneID()
         // ].whichFace(masterPatchStart + i)] =
       //        masterFaceNormals[i];
       //    }

       // //- exchange parallel data
       // // sum because each face is only on one proc
       // reduce(globalMasterFaceNormals, sumOp<vectorField>());

       // use deformed master normals
       const vectorField& globalMasterFaceNormals =
           masterFaceZonePatch.faceNormals();

       // interpolate master normals to the slave and reverse the direction
       // we can use inverseDistance or GGI for the interpolation
       // GGI is much better is we want the actual master normals
       vectorField globalSlaveFaceNormals
           (slaveFaceZonePatch.size(), vector::zero);
       if (ggiInterpolatorPtr)
     {
       //Info << "interpolating master normals with GGI" << endl;
       globalSlaveFaceNormals =
         -ggiInterpolatorPtr->masterToSlave(globalMasterFaceNormals);
     }
       else // inverse distance
     {
       globalSlaveFaceNormals =
         -masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
         (
          globalMasterFaceNormals
          );
     }
       vectorField actualSlaveFaceNormals =
     slaveFaceZonePatch.faceNormals();

       // get local normals from global
       const label slavePatchStart
     = mesh.boundaryMesh()[slavePatchIndex].start();
       forAll(slaveFaceNormals, facei)
             {
               slaveFaceNormals[facei] =
                 globalSlaveFaceNormals
                 [
                  mesh.faceZones()[slaveFaceZoneID()].whichFace
                  (slavePatchStart + facei)
                  ];

       if (mag(slaveFaceNormals[facei]) < SMALL)
         {
           // interpolation sometimes does not work for far
           // away faces so we will use actual slave normals
           // but doesn't really matter as the face is probably
           // not in contact
           slaveFaceNormals[facei] = actualSlaveFaceNormals[facei];
         }
       else
         {
           // make sure they are unity
           slaveFaceNormals[facei] /= mag(slaveFaceNormals[facei]);
         }
             }
     }
   else
     {
       FatalError << "dirichletNeumann::correct()\n"
          << "nonLinear option " << nonLinear
          << " is unknown"
          << exit(FatalError);
     }
