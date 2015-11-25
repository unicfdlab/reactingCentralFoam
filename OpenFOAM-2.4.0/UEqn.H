    fvVectorMatrix UEqn
    (
	fvm::ddt(rho, U)
	+
	fvm::div(phiPos, U)
	+
	fvm::div(phiNeg, U)
	+
	turbulence->divDevRhoReff(U)
	==
	rho*g
	+
	fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
	solve(UEqn == - gradp);
	fvOptions.correct(U);
	K = 0.5*magSqr(U);
    }
