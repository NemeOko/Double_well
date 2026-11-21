# DoubleWell

**Simulating a chaotic, magnetically tuned piezoelectric energy harvester for ultra-low-frequency vibrations (MATLAB).**

Small wireless sensors and IoT devices often sit on machines, structures or people that vibrate slowly (1 to 20 Hz). A standard piezoelectric cantilever harvester only works well near its resonance, which is usually much higher. This project models a harvester with a magnet pair that cancels most of the beam's stiffness (quasi-zero stiffness, QZS). That creates a double-well potential, so the beam can flip between wells and keep generating voltage at low frequencies.

This is a **theoretical** study from my MEng dissertation at the University of Manchester (2025). A prototype has been designed but not yet tested.

## What the code does
- Builds a distributed-parameter model of a bimorph cantilever (Euler-Bernoulli beam, single-mode expansion, piezoelectric coupling, resistive load).
- Adds the nonlinear magnetic force with a point-dipole approximation.
- Scans the magnet gap to find the one that gives near-zero net stiffness at rest.
- Solves the coupled equations with `ode45` and sweeps the excitation frequency to produce voltage and tip-displacement response curves, plus time histories, phase plots and potential-energy plots.

## Results (from the dissertation)
- Best magnet gap found: **10.4 mm**.
- At ultra-low to low frequencies the QZS harvester produced at least 10x the voltage of an equivalent linear harvester, over a broader band.
- Near the linear harvester's resonance, and above it, the linear design performs as well or better.
- Large-amplitude responses depend on the initial conditions, so forward and reverse sweeps differ (hysteresis).
![Voltage frequency response: QZS (forward and reverse sweeps) vs linear harvester](images/voltage_frf_qzs_vs_linear.png)
 
*Voltage frequency response of the QZS harvester (forward and reverse sweeps) against an equivalent linear harvester. Base excitation: [confirm 0.1 g, 1 g or 2 g from Figure 16 of your report].*
 
![Static potential energy of the beam, magnets and total at the QZS gap](images/static_potential_qzs.png)
 
*Static potential energy at the QZS gap: the magnets flatten the bottom of the well so the beam responds to very small, slow excitation.*
 
## Run it
Requires MATLAB with the Symbolic Math Toolbox. Put the three `.m` files in one folder and run:

```matlab
qzs_harvester_main
```

You will be asked for: magnetic force on/off (1 = on), minimum and maximum frequency (Hz), the frequency step, and the load resistance(s) in ohms. The dissertation results used R = 1000 ohm and base accelerations of 0.1 g, 1 g and 2 g (change `ub2dot` in `piezo_beam_resistor.m`). For a reverse sweep, enter a minimum above the maximum with a negative step.

## Limitations
- Theoretical only. The linear model agrees well with published experiments, but the QZS model was only compared with another group's experiments on a different configuration, and agreement was in trends rather than values.
- Assumes perfectly bonded continuous piezo layers and a point-dipole magnet model, which is less accurate at larger gaps.
- Ignores fatigue and ageing of the piezoelectric material.

## Credits
- Beam model builds on Erturk and Inman (2008, 2009) and Rafique and Bonello (2010).
- Bistable harvester comparison: Stanton, McGehee and Mann (2010), *Physica D* 239(10).
  
## Author
Neme Okoroafor, MEng Mechanical Engineering, University of Manchester.
