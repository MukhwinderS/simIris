# simIris
Simulation tool for the IRIS experiment @ TRIUMF


## Download and Installation ##
	
simIris requires ROOT (https://root.cern.ch).

To install simIris use

	git clone https://github.com/iris-triumf/simIris
	cd simIris
	make


## Running simIris ##

From the main simIris folder, you can run simIris with
 
	./bin/simIris --output=/path/to/to/your/output.file --reaction=/path/to/your/reaction.file --geo=/path/to/your/geometry.file --dedx_dir=/path/to/folder/containing/energy_loss_tables --events=number of events to be simulated [--run=run number]

Configuration files for the reaction and the geometry as well as energy loss tables have to be provided, the output file will be a root file which can be further analysed using physIris.


## The Geometry File ##

The geometry file contains a set of parameters describing the experimental setup of IRIS to be used in the simulation. The following parameters are available:
	
	Bs: Size of beam spot in mm (usually ~2mm)			
	ICPressure: Gas pressure of ionization chamber in Torr (usually 19.5 Torr or lower)
	FOIL: Material of the foil (Ag for all experiments using the H2 or D2 target)			
	TFoil: Thickness of the foil in mg/cm2 (4.816 for current Ag foil)
	AoZFoil: Mass to charge ration of foil (2.3 for Ag foil)
	TARGET: Target material (H or D)		
	TTgt: Target thickness in mg/cm2
	TSD[1-2],TSU,TYY[1-8],TYU[1-8]: Detector thicknesses in um for S3-1, S3-2, upstream S3, and 8 YY1 array segments for the downstream and the upstream array
	DSD,DSU,DYY,DYU: Detector distances from target in mm for S3-1, upstream S3, downstream and upstream YY1 arrays. S3-2 is at a fixed distance to S3-1. Convention: Distances for the upstream detectors have to be negative	
	MASK: Use the target mask? (1=yes, 0=no)
	SHIELD: Use the heat shield? (1=yes, 0=no)
	ORIENTATION: Orientation of the target system (0=foil upstream of target, 1=foil downstream of target)


## The Reaction Parameter File ##

The reaction parameter file contains information about the reaction to be simulated. The following parameters are available:

	N: number of particles in outgoing channel, Reactions with up to 6 particles are supported by simIris
	A,a,B,b,c,d,e,f: Define the particles involved in the reaction: A(a,b)B[+c+d+e+f] 
	R1,W1: Excitation energy and resonance widht of beam-like ejectile in MeV
	R2,W2: Excitation energy and resonance widht of target-like ejectile in MeV
	E: Beam energy in MeV
	SHT: Reaction on the solid H2/D2 target? (1=yes, 0=no) if this is 0 the reaction will take place at middle of silver foil
	foil= foil material (107Ag)  

Each parameter in both the geometry and the reaction parameter file has to be separated from the value by an "="-sign without any spaces, e.g. ICPressure=19.5

## Energy loss tables ##

simIris uses SRIM  energy loss tables to calculate the energy loss of particles going through detector materials or the target. Tables in a usable format can be created using LISE++.

## Output files

The output root file of simIris contains a tree with the following objects:

	yd, yu, sd1, sd2, su, csi: Containers for detector hits. They contain the detector thickness, and information about the simulated detector hit: multiplicity, X-, Y-, Z-position, angls, git segments, and depostited energy.

	tlP, blP: Contain information about the generated particle tracks. Contains angles T(heta) and P(hi) and energies in the lab and center-of-mass system. Also saves the energy loss in the silver foil and target as well as the resulting energy behind the target (FoildE, TrgtdE, Ebt).

	tlPdec1/2, blPdec1/2: Same as blP and tlP, but only used when the heavy particle at the simulated excitation energy is unbound.

	det: Contains the same information as yd, yu, etc. but in the same format as data sorted with treeIris. Simulated data from simIris can therefore analysed using physIris (https://github.com/iris-triumf/physIris) the same way as real experimental data.
