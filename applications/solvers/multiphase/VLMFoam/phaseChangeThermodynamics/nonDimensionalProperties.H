
	// Extract the temperature scalar values (i.e., remove the dimensionality)
	forAll(T, cellI) // Loop over all cells of the internal field (i.e. cells at the boundaries not included)
	{
	  Tdimless[cellI] = T[cellI];
	}

	forAll(mesh.boundaryMesh(), patchID) // Iterate for every patch (patch = boundary here)
	{
	     forAll (mesh.boundaryMesh()[patchID],faceI) // Iterate for every cell on the patch
	     {
		if(mesh.boundaryMesh()[patchID].type() != "empty")
		{
		  	Tdimless.boundaryFieldRef()[patchID][faceI] = T.boundaryField()[patchID][faceI];
		}
	     }
	}

	forAll(rho, cellI) // Loop over all cells of the internal field (i.e. cells at the boundaries not included)
	{
	  rhodimless[cellI] = rho[cellI];
	}

	forAll(mesh.boundaryMesh(), patchID) // Iterate for every patch (patch = boundary here)
	{
	     forAll (mesh.boundaryMesh()[patchID],faceI) // Iterate for every cell on the patch
	     {
		if(mesh.boundaryMesh()[patchID].type() != "empty")
		{
		  	rhodimless.boundaryFieldRef()[patchID][faceI] = rho.boundaryField()[patchID][faceI];
		}
	     }
	}
