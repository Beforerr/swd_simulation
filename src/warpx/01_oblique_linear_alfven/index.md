# Oblique linearly polarized Alfvén wave

References:

- @vasquezMakingAlfvenicFluctuation1996

> We now consider hybrid simulations (see Vasquez [1995] for details of the code) of the evolution of an initially linearly polarized . Alfven wave train with A = 0.3 and traveling at OB = 60°. We choose a small wavenumher of 0.025(c/w¢) -1 where c/wm is an ion inertial length. At this wavenumber ion wave dispersion is negligible over the total simulation time. In the simulation box, 4 cycles of this wave will be included. Time in the simulation is normalized to an ion gyrocycle ft. In units of (c/c~)f~, the Alfvdn speed is 1, and the speed of the wave train is cos 60° = 0.5. The wave period is 502 ft-1. Electrons are treated as an isothermal fluid with temperature Te. Ions are assigned a Maxwellian velocity distribution with the bulk speed of an Alfvdn wave and temperature Ti equal to 2",. The total plasma /3 ks taken to be 0.275. We are especially concerned here with small/3 since this condition characterizes re~ons near the Sun. The isothermal sound speed corresponding to this/3 is 0.37, and this value will be used for c, to evaluate Equations (8)-(9). Evaluated in this way, the correspondence with the simulations is close but not exact because ions in a hybrid simulation do not obey a predetermined equation of state.

Notes:

- There seems to have an initialization problem for 1D simulation.

- Polarization is dertemined by the wave dispersion relation. 

    - [ ] Verify the linear polarization

- Numerical heating may be a problem for long time simulation.

    - See issues [Hybrid Code Crash · Issue #4883 · ECP-WarpX/WarpX](https://github.com/ECP-WarpX/WarpX/issues/4883)

    > Hyper resistivity is added to the generalized Ohm's law solver to suppress spurious whistler noise in low density regions.

    - Low theta is unstable

TODO:

- $\gamma$ in `warpx_momentum_spread_expressions`


- More oblique, less anisotropic, more apparent RD-like structure