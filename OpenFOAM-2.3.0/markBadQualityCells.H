{
    Info << "Marking cells with bad quality to limit pressure gradient" << endl;
    
    badQualityCells.resize(0);
    
    const scalar nonOrthogonalityThreshold
    (
	readScalar
	(
	    mesh.solutionDict().subDict("PIMPLE").lookup("nonOrthogonalityThreshold")
	)
    );
    
    const scalar skewnessThreshold
    (
	readScalar
	(
	    mesh.solutionDict().subDict("PIMPLE").lookup("skewnessThreshold")
	)
    );
    
    cellQuality cq (mesh);
    
    scalarField cellNonOrtho (cq.nonOrthogonality());
    
    scalarField cellSkewness (cq.skewness());
    
    forAll(cellNonOrtho, iCell)
    {
	if ( (cellNonOrtho[iCell] >= nonOrthogonalityThreshold) ||
		(cellSkewness[iCell] >= skewnessThreshold) )
	{
	    badQualityCells.append(iCell);
	}
    }
    
    label nTotalBadCells = badQualityCells.size();
    reduce(nTotalBadCells, sumOp<label>());
    
    Info << badQualityCells.size() << " cells marked as bad" << endl
    << "Standard fvm pressure gradients will be applied in this cells instead of KT" << endl;
}
