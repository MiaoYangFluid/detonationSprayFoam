// --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi) == parcels.Srho(rho));
