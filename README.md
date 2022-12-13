# Environment for simulating the AWAKE Run 2 Electron Transfer Line

For a full description of this code and the results see [here](https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.25.101602).

Code for designing, optimising and simulating a 150 MeV electron beamline. Simulation includes real-world
effects such as beam scattering within materials, magnet imperfections, element/beam misalignments, finite detector
resolutions, realistic beamline set-up procedures. 

The simulation code MAD-X was used to model the beam transport, with the bunch tracking simulated using PTC (polymorphic tracking code).
The scattering foil was modelled using Multiple Coulomb Scattering in Python.
Errors were added using Monte-Carlo methods, sampling from Gaussian distributions. 
Encoded constraints include the tunnel width, the placement of plasma cells, the alignment of the laser lines. 

There is functionality to perform both single and multiobjective design optimisation.
Use the flag [pyMOO](https://pymoo.org/) to use genetic algorithms. 

Using **Thin-track** (vs. PTC track) is quicker but less accurate. Similarly, to increase speed track fewer macro-particles, although, again at the expense of accuracy.

### Functionality for optimisation, plotting, error studies, steering, etc. 
For the pseudo-code for these algorithms see Appendix of this [paper](https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.25.101602).

* Top-level script: `runOpt.py`.
* Optimisation environment: `optEnv.py`.
* Multi-objective optimisation environment: `MO_env.py`.
* Error study framework: `errorOptEnv.py`.
* Beam tracking functionality: `get_beam_size.py`.
* Plotting scripts: `plot_save_output.py`.
* Add static errors: `add_errors.madx`.
* Add dynamic errors: `add_errors_jitt.madx`.
* Beamline sequence and model under `general_tt43_python.madx`.


## Optimiser:
* Enter the details of the parameters you want to optimise into x (in runOpt)
* The objective function is under step in OptEnv
* Most of the optimiser parameters are set under runOpt e.g number of particles to track, solver to use, width of foil at beam waist

## ErrorEnv allows for alignment and steering simulations:
- Adds current corrector values to new ones, so can string steering/alignment methods together sequentially
* SVD
* DFS
* One-to-one
* Quad shunting
	* Vary quad strength by 20% and measure at downstream BPM to estimate quad-beam offset
	* Correct for this offset, and iterate along beamline
	* Enter gain and number of iterations of quad shunting to perform
	* Max correction of 100 um
	* Switches higher order magnets on and off (self.switch_on_higher_order())
* Sextupole shunting
* Misc.: 
	* zeroCorr - zero all correctors
	* read_bpms

## Plotting
* Phase vs s (plot.phase())
* Beam distributions (plotmat)
* Heatmap of beam size vs s for varying parameters (heatmap)
* Twiss for chromatic/detuning with amplitude effects (diffTwiss)
* Bar chart of different orders of contributions to beam size (MapTable)
* Layout of beamline (Survey)
* Modelling the plasma cell as quads (plasma_focus)