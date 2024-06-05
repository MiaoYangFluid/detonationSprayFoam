 // --- Solve momentum
        solve
        (
        fvm::ddt(rhoU) + fvc::div(phiUp)
        );

        U.ref() =
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
              ==
               parcels.SU(U)
            );
            rhoU = rho*U;
        }
