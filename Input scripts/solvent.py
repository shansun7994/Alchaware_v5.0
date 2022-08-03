#!/usr/bin/env python
# coding: utf-8
import math
import openmmtools
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
from sys import stdout
import numpy as np
from importlib import reload
from simtk.openmm import Vec3
from simtk import unit
from openmmtools.multistate import MultiStateReporter, MultiStateSampler, ReplicaExchangeSampler, ParallelTemperingSampler, SAMSSampler, ReplicaExchangeAnalyzer
from openmmtools.states import GlobalParameterState
import os

kB = 0.008314472471220214 * unit.kilojoules_per_mole/unit.kelvin
temp = 300*unit.kelvin
kT = kB * temp
kcal = 4.1868 * unit.kilojoules_per_mole
kTtokcal = kT/kcal * unit.kilocalories_per_mole
runflag='run'

#New class to allow interpolation of AB atoms
class ABComposableState(GlobalParameterState):
	lambda_sterics_AB=GlobalParameterState.GlobalParameter('lambda_sterics_AB', standard_value=0.0)
	lambda_electrostatics_AB=GlobalParameterState.GlobalParameter('lambda_electrostatics_AB', standard_value=0.0)

args = sys.argv[1:]
if len(args) > 1:
	raise ValueError('Only take one argument runflag')

if len(args) == 0:
	runflag='run'
else:
	runflag=args[0]
	allowedrunflags = ['run', 'extend', 'recover']
	if runflag not in allowedrunflags:
		raise ValueError('Please select runflag from {}'.format(allowedrunflags))

#Forcefield
proteinXML='/home/ssun/Parameters/protein.ff14SB.xml'
forcefieldXML='/home/ssun/Parameters/gaff-2.11.xml'
waterModel='tip4pew'

#System
ionicStrength=0.15*unit.molar
boxPadding=0.9*unit.nanometers

#Pathway
stericsSteps=5
elecSteps=5
stericsExponent=2

#Dynamics
timeStep=4.0*unit.femtoseconds
stepsPerIteration=1250
productionIterations=1000
equilibrationIterations=50
iterationsPerCheckpoint=100
extendIterations=1000

waterXML='amber14/'+waterModel+'.xml'

print("Padding:", boxPadding)
print("Forcefield:", proteinXML, " ", forcefieldXML, " ", waterXML)

pdb = app.PDBFile('input/solvent.pdb')
forcefield = app.ForceField(proteinXML, waterXML, forcefieldXML, 'input/dual.xml')
solvated = app.Modeller (pdb.topology, pdb.positions)
solvated.addExtraParticles(forcefield)
solvated.addSolvent(forcefield, model=waterModel, ionicStrength=ionicStrength, neutralize=True, padding=boxPadding)
system = forcefield.createSystem(solvated.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005, hydrogenMass=4*unit.amus)
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))
forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
nforces = system.getNumForces()
ligand_a_atoms = [35,36,37,38]
ligand_b_atoms = [39,40,41,42,43,44,45,46,47]
all_ligand_atoms = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
ligand_a_bonds = []
ligand_b_bonds = []
ligand_a_angles =[]
ligand_b_angles =[]
ligand_a_torsions =[]
ligand_b_torsions =[]
ligand_a_types = ['null'] * 48
ligand_b_types = ['null'] * 48
ligand_a_types[0] = 'c3'
ligand_b_types[0] = 'c3'
ligand_a_types[1] = 'n'
ligand_b_types[1] = 'n'
ligand_a_types[2] = 'cz'
ligand_b_types[2] = 'cz'
ligand_a_types[3] = 'nv'
ligand_b_types[3] = 'nv'
ligand_a_types[4] = 'nu'
ligand_b_types[4] = 'nu'
ligand_a_types[5] = 'c3'
ligand_b_types[5] = 'c3'
ligand_a_types[6] = 'ca'
ligand_b_types[6] = 'ca'
ligand_a_types[7] = 'ca'
ligand_b_types[7] = 'ca'
ligand_a_types[8] = 'ca'
ligand_b_types[8] = 'ca'
ligand_a_types[9] = 'ca'
ligand_b_types[9] = 'ca'
ligand_a_types[10] = 'cp'
ligand_b_types[10] = 'cp'
ligand_a_types[11] = 'cp'
ligand_b_types[11] = 'cp'
ligand_a_types[12] = 'ca'
ligand_b_types[12] = 'ca'
ligand_a_types[13] = 'ca'
ligand_b_types[13] = 'ca'
ligand_a_types[14] = 'ca'
ligand_b_types[14] = 'ca'
ligand_a_types[15] = 'ca'
ligand_b_types[15] = 'ca'
ligand_a_types[16] = 'cl'
ligand_b_types[16] = 'cl'
ligand_a_types[17] = 'ca'
ligand_b_types[17] = 'ca'
ligand_a_types[18] = 'ca'
ligand_b_types[18] = 'ca'
ligand_a_types[19] = 'c'
ligand_b_types[19] = 'c'
ligand_a_types[20] = 'o'
ligand_b_types[20] = 'o'
ligand_a_types[21] = 'h1'
ligand_b_types[21] = 'h1'
ligand_a_types[22] = 'h1'
ligand_b_types[22] = 'h1'
ligand_a_types[23] = 'h1'
ligand_b_types[23] = 'h1'
ligand_a_types[24] = 'hn'
ligand_b_types[24] = 'hn'
ligand_a_types[25] = 'hn'
ligand_b_types[25] = 'hn'
ligand_a_types[26] = 'ha'
ligand_b_types[26] = 'ha'
ligand_a_types[27] = 'ha'
ligand_b_types[27] = 'ha'
ligand_a_types[28] = 'ha'
ligand_b_types[28] = 'ha'
ligand_a_types[29] = 'ha'
ligand_b_types[29] = 'ha'
ligand_a_types[30] = 'ha'
ligand_b_types[30] = 'ha'
ligand_a_types[31] = 'ha'
ligand_b_types[31] = 'ha'
ligand_a_types[32] = 'ha'
ligand_b_types[32] = 'ha'
ligand_a_types[33] = 'ha'
ligand_b_types[33] = 'ha'
ligand_a_types[34] = 'hn'
ligand_b_types[34] = 'hn'
ligand_a_types[35] = 'c3'
ligand_b_types[35] = 'c3'
ligand_a_types[36] = 'hc'
ligand_b_types[36] = 'hc'
ligand_a_types[37] = 'hc'
ligand_b_types[37] = 'hc'
ligand_a_types[38] = 'hc'
ligand_b_types[38] = 'hc'
ligand_a_types[39] = 'ca'
ligand_b_types[39] = 'ca'
ligand_a_types[40] = 'ca'
ligand_b_types[40] = 'ca'
ligand_a_types[41] = 'nb'
ligand_b_types[41] = 'nb'
ligand_a_types[42] = 'ca'
ligand_b_types[42] = 'ca'
ligand_a_types[43] = 'nb'
ligand_b_types[43] = 'nb'
ligand_a_types[44] = 'ca'
ligand_b_types[44] = 'ca'
ligand_a_types[45] = 'h4'
ligand_b_types[45] = 'h4'
ligand_a_types[46] = 'h5'
ligand_b_types[46] = 'h5'
ligand_a_types[47] = 'h4'
ligand_b_types[47] = 'h4'
ligand_a_charges = ['null'] * 48
ligand_a_charges[0] = 0.0563
ligand_a_charges[1] = -0.4004
ligand_a_charges[2] = 0.4336
ligand_a_charges[3] = -0.4632
ligand_a_charges[4] = -0.3612
ligand_a_charges[5] = 0.052
ligand_a_charges[6] = -0.1533
ligand_a_charges[7] = -0.122
ligand_a_charges[8] = -0.113
ligand_a_charges[9] = -0.085
ligand_a_charges[10] = -0.007
ligand_a_charges[11] = -0.058
ligand_a_charges[12] = -0.108
ligand_a_charges[13] = -0.113
ligand_a_charges[14] = -0.102
ligand_a_charges[15] = 0.0184
ligand_a_charges[16] = -0.0774
ligand_a_charges[17] = -0.115
ligand_a_charges[18] = -0.1
ligand_a_charges[19] = 0.7005
ligand_a_charges[20] = -0.4675
ligand_a_charges[21] = 0.0847
ligand_a_charges[22] = 0.0847
ligand_a_charges[23] = 0.0847
ligand_a_charges[24] = 0.3327
ligand_a_charges[25] = 0.3297
ligand_a_charges[26] = 0.134
ligand_a_charges[27] = 0.157
ligand_a_charges[28] = 0.161
ligand_a_charges[29] = 0.143
ligand_a_charges[30] = 0.152
ligand_a_charges[31] = 0.161
ligand_a_charges[32] = 0.137
ligand_a_charges[33] = 0.151
ligand_a_charges[34] = 0.3327
ligand_a_charges[35] = -0.0981
ligand_a_charges[36] = 0.080367
ligand_a_charges[37] = 0.080367
ligand_a_charges[38] = 0.080367
ligand_a_charges[39] = -0.3839
ligand_a_charges[40] = 0.4642
ligand_a_charges[41] = -0.671
ligand_a_charges[42] = 0.6679
ligand_a_charges[43] = -0.671
ligand_a_charges[44] = 0.4642
ligand_a_charges[45] = 0.0431
ligand_a_charges[46] = 0.0821
ligand_a_charges[47] = 0.0431
ligand_b_charges = ['null'] * 48
ligand_b_charges[0] = 0.0553
ligand_b_charges[1] = -0.3974
ligand_b_charges[2] = 0.4376
ligand_b_charges[3] = -0.4612
ligand_b_charges[4] = -0.3652
ligand_b_charges[5] = 0.1323
ligand_b_charges[6] = -0.1593
ligand_b_charges[7] = -0.109
ligand_b_charges[8] = -0.114
ligand_b_charges[9] = -0.082
ligand_b_charges[10] = -0.009
ligand_b_charges[11] = -0.058
ligand_b_charges[12] = -0.107
ligand_b_charges[13] = -0.112
ligand_b_charges[14] = -0.101
ligand_b_charges[15] = 0.0174
ligand_b_charges[16] = -0.0774
ligand_b_charges[17] = -0.117
ligand_b_charges[18] = -0.102
ligand_b_charges[19] = 0.7035
ligand_b_charges[20] = -0.4695
ligand_b_charges[21] = 0.086367
ligand_b_charges[22] = 0.086367
ligand_b_charges[23] = 0.086367
ligand_b_charges[24] = 0.3337
ligand_b_charges[25] = 0.3327
ligand_b_charges[26] = 0.141
ligand_b_charges[27] = 0.159
ligand_b_charges[28] = 0.161
ligand_b_charges[29] = 0.144
ligand_b_charges[30] = 0.153
ligand_b_charges[31] = 0.162
ligand_b_charges[32] = 0.135
ligand_b_charges[33] = 0.145
ligand_b_charges[34] = 0.3337
ligand_b_charges[35] = -0.0981
ligand_b_charges[36] = 0.080367
ligand_b_charges[37] = 0.080367
ligand_b_charges[38] = 0.080367
ligand_b_charges[39] = -0.3839
ligand_b_charges[40] = 0.4642
ligand_b_charges[41] = -0.671
ligand_b_charges[42] = 0.6679
ligand_b_charges[43] = -0.671
ligand_b_charges[44] = 0.4642
ligand_b_charges[45] = 0.0431
ligand_b_charges[46] = 0.0821
ligand_b_charges[47] = 0.0431
#Get date from forcefield
bond_force=forces['HarmonicBondForce']
angle_force = forces['HarmonicAngleForce']
torsion_force=forces['PeriodicTorsionForce']
bond_data = forcefield.getGenerators()[0]
angle_data = forcefield.getGenerators()[1]
torsion_data = forcefield.getGenerators()[2]
nonbonded_data = forcefield.getGenerators()[3]
num_nonbonded=len(nonbonded_data.params.paramsForType)
num_bonds=len(bond_data.length)
num_angles=len(angle_data.angle)
num_propers=len(torsion_data.proper)
num_impropers=len(torsion_data.improper)
#Check for abberant bonds that can be caused by an unknown error with water molecules
for index in range(bond_force.getNumBonds()):
        [atom_i, atom_j, r0, K] = bond_force.getBondParameters(index)
        if set([atom_i, atom_j]).intersection(all_ligand_atoms):
                if set([atom_i]).intersection(all_ligand_atoms) and not set([atom_j]).intersection(all_ligand_atoms):
                        print("Fatal Error: bond between ligand",atom_i,"and non ligand",atom_j,". Exiting.")
                        exit(0)
                if set([atom_j]).intersection(all_ligand_atoms) and not set([atom_i]).intersection(all_ligand_atoms):
                        print("Fatal Error: bond between ligand",atom_j,"and non ligand",atom_i,". Exiting.")
                        exit(0)
#Some of the parameters are from ligand A atom types
#Here we correct all these parameters back to those from ligand B atom types
for force_index in range(nforces):
	force = system.getForce(force_index)
	if isinstance(force, mm.NonbondedForce):
		for index in range(force.getNumParticles()):
			if set([index]).intersection(all_ligand_atoms):
				[charge, sigma, epsilon]=force.getParticleParameters(index)
				type=ligand_a_types[index]
				new_charge=ligand_a_charges[index]*unit.elementary_charge
				new_data=nonbonded_data.params.paramsForType.get(type, "none")
				new_sigma=new_data.get("sigma", "none")*unit.nanometers
				new_epsilon=new_data.get("epsilon", "none")*unit.kilojoule_per_mole
				force.setParticleParameters(index, new_charge, new_sigma, new_epsilon)

for force_index in range(nforces):
	bforce = system.getForce(force_index)
	if isinstance(bforce, mm.HarmonicBondForce):
		for index in range(bforce.getNumBonds()):
			[iatom, jatom, r0, K] = bforce.getBondParameters(index)
			if set([iatom,jatom]).intersection(all_ligand_atoms):
				typei=ligand_a_types[iatom]
				typej=ligand_a_types[jatom]
				if set([iatom]).intersection(ligand_b_atoms) or set([jatom]).intersection(ligand_b_atoms):
					typei=ligand_b_types[iatom]
					typej=ligand_b_types[jatom]
				for a in range(num_bonds):
					new_length=bond_data.length[a]*unit.nanometers
					new_force=bond_data.k[a]*unit.kilojoule_per_mole/(unit.nanometers*unit.nanometers)
					type1=bond_data.types1[a]
					type2=bond_data.types2[a]
					if(typei == type1[0] and typej == type2[0] or typei == type2[0] and typej == type1[0]):
						bforce.setBondParameters(index, iatom, jatom, new_length, new_force)
						break
for constraintindex in range(system.getNumConstraints()):
        [iatom, jatom, r0] = system.getConstraintParameters(constraintindex)
        if set([iatom,jatom]).intersection(all_ligand_atoms):
                typei=ligand_a_types[iatom]
                typej=ligand_a_types[jatom]
                if set([iatom]).intersection(ligand_b_atoms) or set([jatom]).intersection(ligand_b_atoms):
                        typei=ligand_b_types[iatom]
                        typej=ligand_b_types[jatom]
                for a in range(num_bonds):
                        new_length=bond_data.length[a]*unit.nanometers
                        type1=bond_data.types1[a]
                        type2=bond_data.types2[a]
                        if(typei == type1[0] and typej == type2[0] or typei == type2[0] and typej == type1[0]):
                                system.setConstraintParameters(constraintindex, iatom, jatom, new_length)
                                break
for force_index in range(nforces):
	aforce = system.getForce(force_index)
	if isinstance(aforce, mm.HarmonicAngleForce):
		for index in range(aforce.getNumAngles()):
			[atom_i, atom_j, atom_k, angle, K] = aforce.getAngleParameters(index)
			if set([atom_i, atom_j, atom_k]).intersection(all_ligand_atoms):
				typei=ligand_a_types[atom_i]
				typej=ligand_a_types[atom_j]
				typek=ligand_a_types[atom_k]
				if set([atom_i]).intersection(ligand_b_atoms) or set([atom_j]).intersection(ligand_b_atoms) or set([atom_k]).intersection(ligand_b_atoms):
					typei=ligand_b_types[atom_i]
					typej=ligand_b_types[atom_j]
					typek=ligand_b_types[atom_k]
				for a in range(num_angles):
					new_angle=angle_data.angle[a]*unit.radian
					new_force=angle_data.k[a]*unit.kilojoule_per_mole/(unit.radian*unit.radian)
					type1=angle_data.types1[a]
					type2=angle_data.types2[a]
					type3=angle_data.types3[a]
					if(typei == type1[0] and typej == type2[0] and typek == type3[0] or typei == type3[0] and typej == type2[0] and typek == type1[0]):
						aforce.setAngleParameters(index, atom_i, atom_j, atom_k, new_angle, new_force)
						break

for force_index in range(nforces):
	tforce = system.getForce(force_index)
	if isinstance(tforce, mm.PeriodicTorsionForce):
		for index in range(tforce.getNumTorsions()):
			[atom_i, atom_j, atom_k, atom_l, periodicity, phase, K] = tforce.getTorsionParameters(index)
			if set([atom_i, atom_j, atom_k, atom_l]).intersection(all_ligand_atoms):
				improper=False
				bonds = list(pdb.topology.bonds())
				for bond in bonds:
					if(atom_i == bond.atom1.index and atom_l == bond.atom2.index or atom_i == bond.atom2.index and atom_l == bond.atom1.index or atom_i == bond.atom1.index and atom_k == bond.atom2.index or atom_i == bond.atom2.index and atom_k == bond.atom1.index or atom_j == bond.atom1.index and atom_l == bond.atom2.index or atom_j == bond.atom2.index and atom_l == bond.atom1.index):
						improper=True
						break
				typei=ligand_a_types[atom_i]
				typej=ligand_a_types[atom_j]
				typek=ligand_a_types[atom_k]
				typel=ligand_a_types[atom_l]
				if set([atom_i]).intersection(ligand_b_atoms) or set([atom_j]).intersection(ligand_b_atoms) or set([atom_k]).intersection(ligand_b_atoms) or set([atom_l]).intersection(ligand_b_atoms):
					typei=ligand_b_types[atom_i]
					typej=ligand_b_types[atom_j]
					typek=ligand_b_types[atom_k]
					typel=ligand_b_types[atom_l]
				if(improper):
					for a in range(num_impropers):
						thisimproper = torsion_data.improper[a]
						new_periodicity=thisimproper.periodicity[0]
						new_phase=thisimproper.phase[0]*unit.radian
						new_force=thisimproper.k[0]*unit.kilojoule_per_mole
						types1=thisimproper.types1
						types2=thisimproper.types2
						types3=thisimproper.types3
						types4=thisimproper.types4
						if (typei in types1) and ((typej in types2 and typek in types3 and typel in types4) or (typej in types2 and typek in types4 and typel in types3) or (typej in types3 and typek in types2 and typel in types4) or (typej in types3 and typek in types4 and typel in types2) or (typej in types4 and typek in types2 and typel in types3) or (typej in types4 and typek in types3 and typel in types2)):
							tforce.setTorsionParameters(index, atom_i, atom_j, atom_k, atom_l, new_periodicity, new_phase, new_force)

				else:
					for a in range(num_propers):
						thisproper = torsion_data.proper[a]
						new_periodicity=thisproper.periodicity[0]
						new_phase=thisproper.phase[0]*unit.radian
						new_force=thisproper.k[0]*unit.kilojoule_per_mole
						types1=thisproper.types1
						types2=thisproper.types2
						types3=thisproper.types3
						types4=thisproper.types4
						if(typei in types1 and typej in types2 and typek in types3 and typel in types4 or typei in types4 and typej in types3 and typek in types2 and typel in types1):
							#Some torsions have multiple k/p/a so only use the right one
							if(new_periodicity != periodicity):
								continue
							tforce.setTorsionParameters(index, atom_i, atom_j, atom_k, atom_l, new_periodicity, new_phase, new_force)

#Identify bonds that are shared between A and B.
#If the parameter types for the bond change as we go from A to B then we turn off the A bond and turn on the B bond
#This is done by creating a new bond for state B
#Search the forcefield for the new parameters.
for index in range(bond_force.getNumBonds()):
	[atom_i, atom_j, r0, K] = bond_force.getBondParameters(index)
	if not set([atom_i, atom_j]).intersection(ligand_a_atoms):
		if not set([atom_i, atom_j]).intersection(ligand_b_atoms):
			#This bond is made of all AB atoms
			type_i_ligand_a=ligand_a_types[atom_i]
			type_i_ligand_b=ligand_b_types[atom_i]
			type_j_ligand_a=ligand_a_types[atom_j]
			type_j_ligand_b=ligand_b_types[atom_j]
			if(type_i_ligand_a != type_i_ligand_b or type_j_ligand_a != type_j_ligand_b):
				length_ligand_b=0
				K_ligand_b=0
				for a in range(num_bonds):
					length=bond_data.length[a]*unit.nanometers
					force=bond_data.k[a]*unit.kilojoule_per_mole/(unit.nanometers*unit.nanometers)
					type1=bond_data.types1[a]
					type2=bond_data.types2[a]
					if (type1[0] == type_i_ligand_b and type2[0] == type_j_ligand_b) or (type1[0] == type_j_ligand_b and type2[0] == type_i_ligand_b):
						length_ligand_b=length
						K_ligand_b=force
						break
				if(length_ligand_b == 0 or K_ligand_b == 0):
					print("No bond defined. Exiting.")
					exit(0)
				if(length_ligand_b == r0 and K_ligand_b == K):
					#No need to change
					continue
				newindex=bond_force.getNumBonds()
				bond_force.addBond(atom_i, atom_j, length_ligand_b, K_ligand_b)
				ligand_a_bonds.append(index)
				ligand_b_bonds.append(newindex)
#Identify angles that are shared between A and B.
#If the parameter types for the angle change as we go from A to B then we turn off the A angles and turn on the B angles
#This is done by creating a new angle for state B
#Search the forcefield for the new parameters.
for index in range(angle_force.getNumAngles()):
	[atom_i, atom_j, atom_k, angle_ligand_a, K_ligand_a] = angle_force.getAngleParameters(index)
	if not set([atom_i, atom_j, atom_k]).intersection(ligand_a_atoms):
		if not set([atom_i, atom_j, atom_k]).intersection(ligand_b_atoms):
			#This angle is made of all AB atoms
			type_i_ligand_a=ligand_a_types[atom_i]
			type_i_ligand_b=ligand_b_types[atom_i]
			type_j_ligand_a=ligand_a_types[atom_j]
			type_j_ligand_b=ligand_b_types[atom_j]
			type_k_ligand_a=ligand_a_types[atom_k]
			type_k_ligand_b=ligand_b_types[atom_k]
			if(type_i_ligand_a != type_i_ligand_b or type_j_ligand_a != type_j_ligand_b or type_k_ligand_a != type_k_ligand_b):
				angle_ligand_b=0
				K_ligand_b=0
				for a in range(num_angles):
					angle=angle_data.angle[a]
					force=angle_data.k[a]
					type1=angle_data.types1[a]
					type2=angle_data.types2[a]
					type3=angle_data.types3[a]
					if (type1[0] == type_i_ligand_b and type2[0] == type_j_ligand_b and type3[0] == type_k_ligand_b) or (type1[0] == type_k_ligand_b and type2[0] == type_j_ligand_b and type3[0] == type_i_ligand_b):
						angle_ligand_b=angle*unit.radian
						K_ligand_b=force*unit.kilojoule_per_mole/(unit.radian*unit.radian)
						break
				if(angle_ligand_b == 0 or K_ligand_b == 0):
					print("No angle defined. Exiting.")
					exit(0)
				if(angle_ligand_b == angle_ligand_a and K_ligand_b == K_ligand_a):
					#No need to change
					continue
				newindex=angle_force.getNumAngles()
				angle_force.addAngle(atom_i, atom_j, atom_k, angle_ligand_b, K_ligand_b)
				ligand_a_angles.append(index)
				ligand_b_angles.append(newindex)

nnforces = system.getNumForces()
for force_index in range(nnforces):
	force = system.getForce(force_index)
	if isinstance(force, mm.HarmonicAngleForce):
			for torsion_index in range(force.getNumAngles()):
					particle1, particle2, particle3, angle, k = force.getAngleParameters(torsion_index)
					if set([particle1, particle2, particle3]).intersection(ligand_a_atoms):
							if set([particle1, particle2, particle3]).intersection(ligand_b_atoms):
									force.setAngleParameters(torsion_index, particle1, particle2, particle3, angle, 0)
	if isinstance(force, mm.PeriodicTorsionForce):
			for torsion_index in range(force.getNumTorsions()):
					particle1, particle2, particle3, particle4, periodicity, phase, k = force.getTorsionParameters(torsion_index)
					if set([particle1, particle2, particle3, particle4]).intersection(ligand_a_atoms):
							if set([particle1, particle2, particle3, particle4]).intersection(ligand_b_atoms):
									force.setTorsionParameters(torsion_index, particle1, particle2, particle3, particle4, periodicity, phase, 0)

previous_atom_i=0
previous_atom_j=0
previous_atom_k=0
previous_atom_l=0
for index in range(torsion_force.getNumTorsions()):
	[atom_i, atom_j, atom_k, atom_l, periodicity, phase, K] = torsion_force.getTorsionParameters(index)
	if(atom_i == previous_atom_i and atom_j == previous_atom_j and atom_k == previous_atom_k and atom_l == previous_atom_l):
		#We have already dealt with this torsion
		previous_atom_i=atom_i
		previous_atom_j=atom_j
		previous_atom_k=atom_k
		previous_atom_l=atom_l
		continue

	if not set([atom_j, atom_k]).intersection(ligand_a_atoms):
		if not set([atom_j, atom_k]).intersection(ligand_b_atoms):
			#Atoms 2 and 3 are AB atoms
			if not set([atom_i, atom_l]).intersection(ligand_a_atoms):
				if not set([atom_i, atom_l]).intersection(ligand_b_atoms):
					#This torsion is made of all AB atoms
					type_i_ligand_a=ligand_a_types[atom_i]
					type_i_ligand_b=ligand_b_types[atom_i]
					type_j_ligand_a=ligand_a_types[atom_j]
					type_j_ligand_b=ligand_b_types[atom_j]
					type_k_ligand_a=ligand_a_types[atom_k]
					type_k_ligand_b=ligand_b_types[atom_k]
					type_l_ligand_a=ligand_a_types[atom_l]
					type_l_ligand_b=ligand_b_types[atom_l]
					if(type_i_ligand_a != type_i_ligand_b or type_j_ligand_a != type_j_ligand_b or type_k_ligand_a != type_k_ligand_b or type_l_ligand_a != type_l_ligand_b):
						periodicity_ligand_b=0
						phase_ligand_b=0
						K_ligand_b=0
						improper=False
						bonds = list(pdb.topology.bonds())
						for bond in bonds:
							if(atom_i == bond.atom1.index and atom_l == bond.atom2.index or atom_i == bond.atom2.index and atom_l == bond.atom1.index or atom_i == bond.atom1.index and atom_k == bond.atom2.index or atom_i == bond.atom2.index and atom_k == bond.atom1.index or atom_j == bond.atom1.index and atom_l == bond.atom2.index or atom_j == bond.atom2.index and atom_l == bond.atom1.index):
								improper=True
								break

						existing_torsion_index=-1
						new_torsion_index=-1
						if(improper):
							previous_atom_i=atom_i
							previous_atom_j=atom_j
							previous_atom_k=atom_k
							previous_atom_l=atom_l
							continue
						else:
							for a in range(num_propers):
								thisproper = torsion_data.proper[a]
								types1 = thisproper.types1
								types2 = thisproper.types2
								types3 = thisproper.types3
								types4 = thisproper.types4
								if (type_i_ligand_a in types1 and type_j_ligand_a in types2 and type_k_ligand_a in types3 and type_l_ligand_a in types4) or (type_i_ligand_a in types4 and type_j_ligand_a in types3 and type_k_ligand_a in types2 and type_l_ligand_a in types1):
									existing_torsion_index=a
								if (type_i_ligand_b in types1 and type_j_ligand_b in types2 and type_k_ligand_b in types3 and type_l_ligand_b in types4) or (type_i_ligand_b in types4 and type_j_ligand_b in types3 and type_k_ligand_b in types2 and type_l_ligand_b in types1):
									new_torsion_index=a


						if(existing_torsion_index == -1 or new_torsion_index == -1):
							print(existing_torsion_index, new_torsion_index)
							print(type_i_ligand_a, type_j_ligand_a, type_k_ligand_a, type_l_ligand_a)
							print(type_i_ligand_b, type_j_ligand_b, type_k_ligand_b, type_l_ligand_b)
							print("No torsion parameters found. Perhaps you forgot to add the default torsion in the GAFF file? Exiting.")
							exit(0)

						existingtorsion=torsion_data.proper[existing_torsion_index]
						newtorsion=torsion_data.proper[new_torsion_index]

						for e in range(len(existingtorsion.k)):
							ligand_a_torsions.append(index+e)
						for n in range(len(newtorsion.k)):
							newindex=torsion_force.getNumTorsions()
							torsion_force.addTorsion(atom_i, atom_j, atom_k, atom_l, newtorsion.periodicity[n], newtorsion.phase[n]*unit.radian, newtorsion.k[n]*unit.kilojoule_per_mole)
							ligand_b_torsions.append(newindex)
	previous_atom_i=atom_i
	previous_atom_j=atom_j
	previous_atom_k=atom_k
	previous_atom_l=atom_l

#Keep one AB-AB-A-A torsion and one AB-AB-B-B torsion. Scale all other AB-AB-A-A and AB-AB-B-B torsions
#Keep all AB-A-A-A and AB-B-B-B torsions
#Scale all AB-AB-AB-A and AB-AB-AB-B torsions
fixed_single_a_torsion=False
fixed_single_b_torsion=False
indices=iter(range(0,torsion_force.getNumTorsions()))
for index in indices:
	[atom_i, atom_j, atom_k, atom_l, periodicity, phase, K] = torsion_force.getTorsionParameters(index)
	improper=False
	bonds = list(pdb.topology.bonds())
	for bond in bonds:
		if(atom_i == bond.atom1.index and atom_l == bond.atom2.index or atom_i == bond.atom2.index and atom_l == bond.atom1.index or atom_i == bond.atom1.index and atom_k == bond.atom2.index or atom_i == bond.atom2.index and atom_k == bond.atom1.index or atom_j == bond.atom1.index and atom_l == bond.atom2.index or atom_j == bond.atom2.index and atom_l == bond.atom1.index):
			improper=True
			break
	if(improper):
		continue
	else:
		if set([atom_i]).intersection(ligand_a_atoms):
			if set([atom_j]).intersection(ligand_a_atoms):
				if set([atom_k]).intersection(ligand_a_atoms):
					if not set([atom_l]).intersection(ligand_a_atoms):
						continue
				else:
					if not set([atom_l]).intersection(ligand_b_atoms):
						if(fixed_single_a_torsion == False):
							fixed_single_a_torsion=True
							print("Special Prop A - ", atom_i, atom_j, atom_k, atom_l)
							if(index+1 < torsion_force.getNumTorsions()):
								[atom_i2, atom_j2, atom_k2, atom_l2, periodicity, phase, K] = torsion_force.getTorsionParameters(index+1)
								if(atom_i == atom_i2 and atom_j == atom_j2 and atom_k == atom_k2 and atom_l == atom_l2):
									next(indices, None)
									if(index+2 < torsion_force.getNumTorsions()):
										[atom_i3, atom_j3, atom_k3, atom_l3, periodicity, phase, K] = torsion_force.getTorsionParameters(index+2)
										if(atom_i==atom_i3 and atom_j==atom_j3 and atom_k==atom_k3 and atom_l==atom_l3):
											next(indices, None)
											if(index+3 < torsion_force.getNumTorsions()):
												[atom_i4, atom_j4, atom_k4, atom_l4, periodicity, phase, K] = torsion_force.getTorsionParameters(index+3)
												if(atom_i==atom_i4 and atom_j==atom_j4 and atom_k==atom_k4 and atom_l==atom_l4):
													next(indices, None)
							continue
						else:
							ligand_a_torsions.append(index)
							print("Prop A - ", atom_i, atom_j, atom_k, atom_l)
			else:
				if not set([atom_k]).intersection(ligand_a_atoms) and not set([atom_k]).intersection(ligand_b_atoms):
					if set([atom_k]).intersection(ligand_a_atoms):
						#Check this one
						ligand_a_torsions.append(index)
						print("Prop A - ", atom_i, atom_j, atom_k, atom_l)
					elif not set([atom_k]).intersection(ligand_b_atoms):
						ligand_a_torsions.append(index)
						print("Prop A - ", atom_i, atom_j, atom_k, atom_l)
		elif set([atom_i]).intersection(ligand_b_atoms):
			if set([atom_j]).intersection(ligand_b_atoms):
				if set([atom_k]).intersection(ligand_b_atoms):
					if not set([atom_l]).intersection(ligand_b_atoms):
						continue
				else:
					if not set([atom_l]).intersection(ligand_a_atoms):
						if(fixed_single_b_torsion == False):
							fixed_single_b_torsion=True
							print("Special Prop B - ", atom_i, atom_j, atom_k, atom_l)
							if(index+1 < torsion_force.getNumTorsions()):
								[atom_i2, atom_j2, atom_k2, atom_l2, periodicity, phase, K] = torsion_force.getTorsionParameters(index+1)
								if(atom_i==atom_i2 and atom_j==atom_j2 and atom_k==atom_k2 and atom_l==atom_l2):
									next(indices, None)
									if(index+2 < torsion_force.getNumTorsions()):
										[atom_i3, atom_j3, atom_k3, atom_l3, periodicity, phase, K] = torsion_force.getTorsionParameters(index+2)
										if(atom_i==atom_i3 and atom_j==atom_j3 and atom_k==atom_k3 and atom_l==atom_l3):
											next(indices, None)
											if(index+3 < torsion_force.getNumTorsions()):
												[atom_i4, atom_j4, atom_k4, atom_l4, periodicity, phase, K] = torsion_force.getTorsionParameters(index+3)
												if(atom_i==atom_i4 and atom_j==atom_j4 and atom_k==atom_k4 and atom_l==atom_l4):
													next(indices, None)
							continue
						else:
							ligand_b_torsions.append(index)
							print("Prop B - ", atom_i, atom_j, atom_k, atom_l)
			else:
				if not set([atom_k]).intersection(ligand_a_atoms) and not set([atom_k]).intersection(ligand_b_atoms):
					if set([atom_k]).intersection(ligand_b_atoms):
						#Check this one
						ligand_b_torsions.append(index)
						print("Prop B - ", atom_i, atom_j, atom_k, atom_l)
					elif not set([atom_k]).intersection(ligand_a_atoms):
						ligand_b_torsions.append(index)
						print("Prop B - ", atom_i, atom_j, atom_k, atom_l)
		else:
			if set([atom_j]).intersection(ligand_a_atoms):
				continue
			elif set([atom_j]).intersection(ligand_b_atoms):
				continue
			else:
				if set([atom_k]).intersection(ligand_a_atoms):
					if(fixed_single_a_torsion == False):
						fixed_single_a_torsion=True
						print("Special Prop A - ", atom_i, atom_j, atom_k, atom_l)
						if(index+1 < torsion_force.getNumTorsions()):
							[atom_i2, atom_j2, atom_k2, atom_l2, periodicity, phase, K] = torsion_force.getTorsionParameters(index+1)
							if(atom_i==atom_i2 and atom_j==atom_j2 and atom_k==atom_k2 and atom_l==atom_l2):
								next(indices, None)
								if(index+2 < torsion_force.getNumTorsions()):
									[atom_i3, atom_j3, atom_k3, atom_l3, periodicity, phase, K] = torsion_force.getTorsionParameters(index+2)
									if(atom_i==atom_i3 and atom_j==atom_j3 and atom_k==atom_k3 and atom_l==atom_l3):
										next(indices, None)
										if(index+3 < torsion_force.getNumTorsions()):
											[atom_i4, atom_j4, atom_k4, atom_l4, periodicity, phase, K] = torsion_force.getTorsionParameters(index+3)
											if(atom_i==atom_i4 and atom_j==atom_j4 and atom_k==atom_k4 and atom_l==atom_l4):
												next(indices, None)
						continue
					else:
						ligand_a_torsions.append(index)
						print("Prop A - ", atom_i, atom_j, atom_k, atom_l)
				elif set([atom_k]).intersection(ligand_b_atoms):
					if(fixed_single_b_torsion == False):
						fixed_single_b_torsion=True
						print("Special Prop B - ", atom_i, atom_j, atom_k, atom_l)
						if(index+1 < torsion_force.getNumTorsions()):
							[atom_i2, atom_j2, atom_k2, atom_l2, periodicity, phase, K] = torsion_force.getTorsionParameters(index+1)
							if(atom_i==atom_i2 and atom_j==atom_j2 and atom_k==atom_k2 and atom_l==atom_l2):
								next(indices, None)
								if(index+2 < torsion_force.getNumTorsions()):
									[atom_i3, atom_j3, atom_k3, atom_l3, periodicity, phase, K] = torsion_force.getTorsionParameters(index+2)
									if(atom_i==atom_i3 and atom_j==atom_j3 and atom_k==atom_k3 and atom_l==atom_l3):
										next(indices, None)
										if(index+3 < torsion_force.getNumTorsions()):
											[atom_i4, atom_j4, atom_k4, atom_l4, periodicity, phase, K] = torsion_force.getTorsionParameters(index+3)
											if(atom_i==atom_i4 and atom_j==atom_j4 and atom_k==atom_k4 and atom_l==atom_l4):
												next(indices, None)
						continue
					else:
						ligand_b_torsions.append(index)
						print("Prop B - ", atom_i, atom_j, atom_k, atom_l)
				else:
					if set([atom_l]).intersection(ligand_a_atoms):
						ligand_a_torsions.append(index)
						print("Prop A - ", atom_i, atom_j, atom_k, atom_l)
					elif set([atom_l]).intersection(ligand_b_atoms):
						ligand_b_torsions.append(index)
						print("Prop B - ", atom_i, atom_j, atom_k, atom_l)

num_a_atoms=len(ligand_a_atoms)
num_b_atoms=len(ligand_b_atoms)
num_a_bonds=len(ligand_a_bonds)
num_b_bonds=len(ligand_b_bonds)
num_a_angles=len(ligand_a_angles)
num_b_angles=len(ligand_b_angles)
num_a_torsions=len(ligand_a_torsions)
num_b_torsions=len(ligand_b_torsions)

#Setup alchemical system
reload(openmmtools.alchemy)
factory = openmmtools.alchemy.AbsoluteAlchemicalFactory(consistent_exceptions=False, split_alchemical_forces = True, alchemical_pme_treatment = 'exact')
reference_system = system

#OpenMM crashes if one adds nonexistent alchemical angles or torsions
if(num_a_bonds == 0):
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms,alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')
else:
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')

if(num_b_bonds == 0):
	if(num_b_angles == 0):
			if(num_b_torsions == 0):
					alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, name='B')
			else:
					alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms,alchemical_torsions = ligand_b_torsions, name='B')
	else:
			if(num_b_torsions == 0):
					alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_angles = ligand_b_angles, name='B')
			else:
					alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_angles = ligand_b_angles, alchemical_torsions = ligand_b_torsions, name='B')
else:
	if(num_b_angles == 0):
		if(num_b_torsions == 0):
			alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_bonds=ligand_b_bonds, name='B')
		else:
			alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_bonds=ligand_b_bonds, alchemical_torsions = ligand_b_torsions, name='B')
	else:
		if(num_b_torsions == 0):
			alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_bonds=ligand_b_bonds, alchemical_angles = ligand_b_angles, name='B')
		else:
			alchemical_region_B = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_b_atoms, alchemical_bonds=ligand_b_bonds, alchemical_angles = ligand_b_angles, alchemical_torsions = ligand_b_torsions, name='B')


alchemical_system_in = factory.create_alchemical_system(reference_system, alchemical_regions = [alchemical_region_A, alchemical_region_B])

#OpenMM creates angles and dihedrals based on bonds.
#This means we have some angles and dihedrals which span A/B which is not correct
#Here we identify and delete these angles and dihedrals
nalchforces = alchemical_system_in.getNumForces()
for force_index in range(nalchforces):
		force = alchemical_system_in.getForce(force_index)
		if isinstance(force, mm.HarmonicAngleForce):
				for torsion_index in range(force.getNumAngles()):
						particle1, particle2, particle3, angle, k = force.getAngleParameters(torsion_index)
						if set([particle1, particle2, particle3]).intersection(ligand_a_atoms):
								if set([particle1, particle2, particle3]).intersection(ligand_b_atoms):
										force.setAngleParameters(torsion_index, particle1, particle2, particle3, angle, 0)
		if isinstance(force, mm.PeriodicTorsionForce):
				for torsion_index in range(force.getNumTorsions()):
						particle1, particle2, particle3, particle4, periodicity, phase, k = force.getTorsionParameters(torsion_index)
						if set([particle1, particle2, particle3, particle4]).intersection(ligand_a_atoms):
								if set([particle1, particle2, particle3, particle4]).intersection(ligand_b_atoms):
										force.setTorsionParameters(torsion_index, particle1, particle2, particle3, particle4, periodicity, phase, 0)

lj14scale=0.5
coulomb14scale=0.8333333333333334
#Fix incorrect 1-4 correction terms
for force_index in range(nalchforces):
        force = alchemical_system_in.getForce(force_index)
        if isinstance(force, mm.CustomBondForce):
                globalparameter=force.getGlobalParameterName(0)
                if(globalparameter != "lambda_sterics_B"):
                        continue
                for index in range(force.getNumBonds()):
                        particle1, particle2, parameters = force.getBondParameters(index)
                        if(set([particle2]).intersection(ligand_b_atoms)):
                                oldsigma1sigma2=parameters[0];
                                oldepsilon1iepsilon2=parameters[1];
                                type=ligand_b_types[particle1]
                                new_data=nonbonded_data.params.paramsForType.get(type, "none")
                                newsigma1=new_data.get("sigma", "none")
                                newepsilon1=new_data.get("epsilon", "none")
                                type=ligand_b_types[particle2]
                                new_data=nonbonded_data.params.paramsForType.get(type, "none")
                                newsigma2=new_data.get("sigma", "none")
                                newepsilon2=new_data.get("epsilon", "none")
                                newsigma1sigma2=0.5*(newsigma1+newsigma2)
                                newepsilon1iepsilon2=lj14scale*math.sqrt(newepsilon1*newepsilon2)
                                newparameters=[newsigma1sigma2,newepsilon1iepsilon2]
                                force.setBondParameters(index, particle1, particle2, newparameters)
        if isinstance(force, mm.NonbondedForce):
                for index in range(force.getNumExceptionParameterOffsets()):
                        parameter, exception, charge, sigma, epsilon = force.getExceptionParameterOffset(index)
                        if(parameter == "lambda_electrostatics_B"):
                                particle1, particle2, echarge, esigma, eepsilon = force.getExceptionParameters(exception)
                                newq1=ligand_b_charges[particle1]*unit.elementary_charge
                                newq2=ligand_b_charges[particle2]*unit.elementary_charge
                                newq1q2=coulomb14scale*newq1*newq2
                                force.setExceptionParameterOffset(index, parameter, exception, newq1q2, sigma, epsilon)

#Create CustomBondForces for any fixes to the 1-4 sterics that are needed. They may remain empty.
region_A_sterics14_energy_expression = ('U_sterics;'
                                                'U_sterics = ((lambda_sterics_A)^softcore_a)*4*epsilon*x*(x-1.0);'
                                                'x = (sigma/reff_sterics)^6;'
                                                'reff_sterics = sigma*((softcore_alpha*(1.0-(lambda_sterics_A))^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);')
region_B_sterics14_energy_expression = ('U_sterics;'
                                                'U_sterics = ((lambda_sterics_B)^softcore_a)*4*epsilon*x*(x-1.0);'
                                                'x = (sigma/reff_sterics)^6;'
                                                'reff_sterics = sigma*((softcore_alpha*(1.0-(lambda_sterics_B))^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);')

region_A_sterics14_force = mm.CustomBondForce(region_A_sterics14_energy_expression)
region_A_sterics14_force.addPerBondParameter("sigma")
region_A_sterics14_force.addPerBondParameter("epsilon")
region_A_sterics14_force.addGlobalParameter('lambda_sterics_A', 1)
region_A_sterics14_force.addGlobalParameter('softcore_alpha', .5)
region_A_sterics14_force.addGlobalParameter('softcore_beta', 0)
region_A_sterics14_force.addGlobalParameter('softcore_a', 1)
region_A_sterics14_force.addGlobalParameter('softcore_b', 1)
region_A_sterics14_force.addGlobalParameter('softcore_c', 6)
region_A_sterics14_force.addGlobalParameter('softcore_d', 1)
region_A_sterics14_force.addGlobalParameter('softcore_e', 1)
region_A_sterics14_force.addGlobalParameter('softcore_f', 2)

region_B_sterics14_force = mm.CustomBondForce(region_B_sterics14_energy_expression)
region_B_sterics14_force.addPerBondParameter("sigma")
region_B_sterics14_force.addPerBondParameter("epsilon")
region_B_sterics14_force.addGlobalParameter('lambda_sterics_B', 1)
region_B_sterics14_force.addGlobalParameter('softcore_alpha', .5)
region_B_sterics14_force.addGlobalParameter('softcore_beta', 0)
region_B_sterics14_force.addGlobalParameter('softcore_a', 1)
region_B_sterics14_force.addGlobalParameter('softcore_b', 1)
region_B_sterics14_force.addGlobalParameter('softcore_c', 6)
region_B_sterics14_force.addGlobalParameter('softcore_d', 1)
region_B_sterics14_force.addGlobalParameter('softcore_e', 1)
region_B_sterics14_force.addGlobalParameter('softcore_f', 2)

for force_index in range(nalchforces):
        force = alchemical_system_in.getForce(force_index)
        if isinstance(force, mm.NonbondedForce):
                for index in range(force.getNumExceptions()):
                        particle1, particle2, charge, sigma, epsilon = force.getExceptionParameters(index)
                        zerocharge=0.0*unit.elementary_charge*unit.elementary_charge
                        zeroepsilon=0.0*unit.kilojoule_per_mole
                        zerosigma=0.0*unit.nanometers
                        onesigma=1.0*unit.nanometers
                        if not(set([particle1, particle2]).intersection(all_ligand_atoms)):
                                #Not in the ligand
                                continue
                        if(sigma == onesigma):
                                #1,2 or 1,3
                                continue
                        if(not set([particle1]).intersection(ligand_a_atoms) and not set([particle1]).intersection(ligand_b_atoms)):
                                if(set([particle2]).intersection(ligand_a_atoms)):
                                        particle1_charge_a=ligand_a_charges[particle1]
                                        particle2_charge_a=ligand_a_charges[particle2]
                                        oldq1q2=coulomb14scale*particle1_charge_a*particle2_charge_a
                                        force.addExceptionParameterOffset('lambda_electrostatics_A',index,oldq1q2,zerosigma,zeroepsilon)
                                        force.setExceptionParameters(index,particle1,particle2,zerocharge,sigma,epsilon)
                                        print("Fixed elec exception AB-A",index,particle1,particle2,oldq1q2)
                                        continue
                                elif(set([particle2]).intersection(ligand_b_atoms)):
                                        particle1_charge_b=ligand_b_charges[particle1]
                                        particle2_charge_b=ligand_b_charges[particle2]
                                        newq1q2=coulomb14scale*particle1_charge_b*particle2_charge_b
                                        force.addExceptionParameterOffset('lambda_electrostatics_B',index,newq1q2,zerosigma,zeroepsilon)
                                        force.setExceptionParameters(index,particle1,particle2,zerocharge,sigma,epsilon)
                                        print("Fixed elec exception AB-B",index,particle1,particle2,newq1q2)
                                        continue

                        particle1_charge_a=ligand_a_charges[particle1]
                        particle1_charge_b=ligand_b_charges[particle1]
                        particle2_charge_a=ligand_a_charges[particle2]
                        particle2_charge_b=ligand_b_charges[particle2]

                        #Alwaysfixcharges
                        oldq1q2=coulomb14scale*particle1_charge_a*particle2_charge_a
                        newq1q2=coulomb14scale*particle1_charge_b*particle2_charge_b

                        if(oldq1q2 != newq1q2):
                            force.addExceptionParameterOffset('lambda_electrostatics_A',index,oldq1q2,zerosigma,zeroepsilon)
                            force.addExceptionParameterOffset('lambda_electrostatics_B',index,newq1q2,zerosigma,zeroepsilon)
                            force.setExceptionParameters(index,particle1,particle2,zerocharge,sigma,epsilon)

                        particle1_type_a=ligand_a_types[particle1]
                        particle1_type_b=ligand_b_types[particle1]
                        particle2_type_a=ligand_a_types[particle2]
                        particle2_type_b=ligand_b_types[particle2]

                        particle1_data_a=nonbonded_data.params.paramsForType.get(particle1_type_a,"none")
                        particle1_data_b=nonbonded_data.params.paramsForType.get(particle1_type_b,"none")
                        particle2_data_a=nonbonded_data.params.paramsForType.get(particle2_type_a,"none")
                        particle2_data_b=nonbonded_data.params.paramsForType.get(particle2_type_b,"none")

                        particle1_epsilon_a=particle1_data_a.get("epsilon","none")
                        particle1_epsilon_b=particle1_data_b.get("epsilon","none")
                        particle2_epsilon_a=particle2_data_a.get("epsilon","none")
                        particle2_epsilon_b=particle2_data_b.get("epsilon","none")

                        particle1_sigma_a=particle1_data_a.get("sigma","none")
                        particle1_sigma_b=particle1_data_b.get("sigma","none")
                        particle2_sigma_a=particle2_data_a.get("sigma","none")
                        particle2_sigma_b=particle2_data_b.get("sigma","none")

                        if(particle1_epsilon_a != particle1_epsilon_b or particle2_epsilon_a != particle2_epsilon_b or particle1_sigma_a != particle1_sigma_b or particle2_sigma_a != particle2_sigma_b):
                                olde1e2=lj14scale*math.sqrt(particle1_epsilon_a*particle2_epsilon_a)
                                newe1e2=lj14scale*math.sqrt(particle1_epsilon_b*particle2_epsilon_b)
                                olds1s2=0.5*(particle1_sigma_a+particle2_sigma_a)
                                news1s2=0.5*(particle1_sigma_b+particle2_sigma_b)

                                region_A_sterics14_force.addBond(particle1, particle2, [olds1s2, olde1e2])
                                region_B_sterics14_force.addBond(particle1, particle2, [news1s2, newe1e2])
                                #If we have created an offset for the charge we must set it to zero
                                if(oldq1q2 != newq1q2):
                                    force.setExceptionParameters(index, particle1, particle2, zerocharge, sigma, zeroepsilon)
                                else:
                                    force.setExceptionParameters(index, particle1, particle2, oldq1q2, sigma, zeroepsilon)

new_force_index=alchemical_system_in.getNumForces()
region_A_sterics14_force.setForceGroup(new_force_index)
alchemical_system_in.addForce(region_A_sterics14_force)
new_force_index=alchemical_system_in.getNumForces()
region_B_sterics14_force.setForceGroup(new_force_index)
alchemical_system_in.addForce(region_B_sterics14_force)
#Create new lambda variables to control the shared atoms.
#We use this to interpolate between the parameters in A and B
alchemical_forces = {alchemical_system_in.getForce(index).__class__.__name__: alchemical_system_in.getForce(index) for index in range(alchemical_system_in.getNumForces())}
nonbonded_force = alchemical_forces['NonbondedForce']
nonbonded_force.addGlobalParameter('lambda_sterics_AB', 1.0)
nonbonded_force.addGlobalParameter('lambda_electrostatics_AB', 1.0)

alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
alchemical_state_B = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'B')
alchemical_state_AB = ABComposableState(lambda_sterics_AB=0, lambda_electrostatics_AB=0)
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=300*unit.kelvin, pressure=1*unit.bar)
composable_states = [alchemical_state_A, alchemical_state_B, alchemical_state_AB]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)
integrator=mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
context = compound_state.create_context(integrator)
alchemical_system_in=context.getSystem()
for force in alchemical_system_in.getForces():
        if isinstance(force, mm.CustomNonbondedForce):
                globalparameter=force.getGlobalParameterName(0)
                if(globalparameter == "lambda_sterics_A"):
                        force.addGlobalParameter('lambda_sterics_AB', 1.0)
                        force.addPerParticleParameter("sigma_offset")
                        force.addPerParticleParameter("epsilon_offset")
                        for particle_index in range(force.getNumParticles()):
                                [a,b] = force.getParticleParameters(particle_index)
                                force.setParticleParameters(particle_index, [a,b,0,0])
                        force.setEnergyFunction("U_sterics;U_sterics = ((lambda_sterics_A)^softcore_a)*4*epsilon*x*(x-1.0);x = (sigma/reff_sterics)^6;reff_sterics = sigma*((softcore_alpha*(1.0-(lambda_sterics_A))^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);epsilon = sqrt((epsilon1+epsilon_offset1)*(epsilon2+epsilon_offset2));sigma = 0.5*((sigma1+sigma_offset1) + (sigma2+sigma_offset2));")
                        force.updateParametersInContext(context)
                        custom_nonbonded_force_A=force
                elif(globalparameter == "lambda_sterics_B"):
                        force.addGlobalParameter('lambda_sterics_AB', 1.0)
                        force.addPerParticleParameter("sigma_offset")
                        force.addPerParticleParameter("epsilon_offset")
                        for particle_index in range(force.getNumParticles()):
                                [a,b] = force.getParticleParameters(particle_index)
                                force.setParticleParameters(particle_index, [a,b,0,0])
                        force.setEnergyFunction("U_sterics;U_sterics = ((lambda_sterics_B)^softcore_a)*4*epsilon*x*(x-1.0);x = (sigma/reff_sterics)^6;reff_sterics = sigma*((softcore_alpha*(1.0-(lambda_sterics_B))^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);epsilon = sqrt((epsilon1+epsilon_offset1)*(epsilon2+epsilon_offset2));sigma = 0.5*((sigma1+sigma_offset1) + (sigma2+sigma_offset2));")
                        force.updateParametersInContext(context)
                        custom_nonbonded_force_B=force

alchemical_forces = {alchemical_system_in.getForce(index).__class__.__name__: alchemical_system_in.getForce(index) for index in range(alchemical_system_in.getNumForces())}
nonbonded_force = alchemical_forces['NonbondedForce']

#Use offsets to interpolate
nonbonded_data = forcefield.getGenerators()[3]
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 0, -0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 1, 0.003, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 2, 0.004, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 3, 0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 4, -0.004, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 5, 0.0803, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 6, -0.006, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 7, 0.013, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 8, -0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 9, 0.003, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 10, -0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 12, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 13, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 14, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 15, -0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 17, -0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 18, -0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 19, 0.003, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 20, -0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 21, 0.001667, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 22, 0.001667, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 23, 0.001667, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 24, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 25, 0.003, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 26, 0.007, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 27, 0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 29, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 30, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 31, 0.001, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 32, -0.002, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 33, -0.006, 0, 0)
nonbonded_force.addParticleParameterOffset('lambda_electrostatics_AB', 34, 0.001, 0, 0)
alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
alchemical_state_B = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'B')
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=300*unit.kelvin, pressure=1*unit.bar)
composable_states = [alchemical_state_A, alchemical_state_B, alchemical_state_AB]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)

stericsSteps+=1
elecSteps+=1
nstates=(2*stericsSteps)+elecSteps-2
print("There will be ", nstates, " states in total")
print("Lambda exponent: ", stericsExponent)
print("stepsPerIteration:", stepsPerIteration, " productionIterations: ", productionIterations, "equilibrationIterations: ", equilibrationIterations)
print("Timestep: ", timeStep)
box_vec = alchemical_system_in.getDefaultPeriodicBoxVectors()
print("Box vectors:", box_vec)

#Setup lambda schedule
lambdas = []
for j in range(nstates):
		column = []
		for i in range(2):
				column.append(1)
		for i in range(4):
				column.append(0)
		lambdas.append(column)

lambdas_sterics_A = np.linspace(1.0, 0.0, stericsSteps)
lambdas_elec_A = np.linspace(1.0, 0.0, elecSteps)
lambdas_sterics_AB = np.linspace(0.0, 0.5, stericsSteps)
lambdas_elec_AB = np.linspace(0.0, 1.0, elecSteps)
lambdas_sterics_B = np.linspace(0.0, 1.0, stericsSteps)
lambdas_elec_B = np.linspace(0.0, 1.0, elecSteps)

for j in range(stericsSteps):
	lambdas_sterics_A[j]=pow(lambdas_sterics_A[j],stericsExponent)
	lambdas_sterics_B[j]=pow(lambdas_sterics_B[j],stericsExponent)

#First switch on B sterics and AB sterics halfway
for j in range(stericsSteps):
	lambdas[j][4]=lambdas_sterics_B[j]
	lambdas[j][2]=lambdas_sterics_AB[j]

#B sterics and AB sterics stay
for j in range(elecSteps+stericsSteps-1):
	lambdas[j+stericsSteps-1][4]=1.0
	lambdas[j+stericsSteps-1][2]=0.5

#Then swap elec
for j in range(elecSteps):
	lambdas[j+stericsSteps-1][1]=lambdas_elec_A[j]
	lambdas[j+stericsSteps-1][3]=lambdas_elec_AB[j]
	lambdas[j+stericsSteps-1][5]=lambdas_elec_B[j]

#They stay swapped
for j in range(stericsSteps):
	lambdas[j+stericsSteps+elecSteps-2][1]=0
	lambdas[j+stericsSteps+elecSteps-2][3]=1
	lambdas[j+stericsSteps+elecSteps-2][5]=1

#Last switch off A torsions and AB sterics from halfway
for j in range(stericsSteps):
	lambdas[j+stericsSteps+elecSteps-2][0]=lambdas_sterics_A[j]
	lambdas[j+stericsSteps+elecSteps-2][2]=0.5+lambdas_sterics_AB[j]

#Sanity check
print("")
print("Lambdas matrix")
print("Lsterics_A      Lelec_A         Lsterics_AB     Lelec_AB        Lsterics_B      Lelec_B")
for j in range(len(lambdas)):
	for i in range(len(lambdas[j])):
		print("%-15.2f" % lambdas[j][i], end=' ')
	print("")
print("")

sampler_states = list()
thermodynamic_states = list()

for k in range(nstates):
    compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
    if(num_a_atoms != 0):
        compound_state.lambda_sterics_A=lambdas[k][0]
        compound_state.lambda_electrostatics_A=lambdas[k][1]
    if(num_a_bonds != 0):
        compound_state.lambda_bonds_A=lambdas[k][1]
    if(num_a_angles != 0):
        compound_state.lambda_angles_A=lambdas[k][1]
    if(num_a_torsions != 0):
        compound_state.lambda_torsions_A=lambdas[k][1]
    compound_state.lambda_sterics_AB=lambdas[k][2]
    compound_state.lambda_electrostatics_AB=lambdas[k][3]
    if(num_b_atoms != 0):
        compound_state.lambda_sterics_B=lambdas[k][4]
        compound_state.lambda_electrostatics_B=lambdas[k][5]
    if(num_b_bonds != 0):
        compound_state.lambda_bonds_B=lambdas[k][5]
    if(num_b_angles != 0):
        compound_state.lambda_angles_B=lambdas[k][5]
    if(num_b_torsions != 0):
        compound_state.lambda_torsions_B=lambdas[k][5]
    sys = compound_state.get_system()
    sampler_states.append(openmmtools.states.SamplerState(positions=solvated.positions, box_vectors=box_vec))
    thermodynamic_states.append(compound_state)

print("Integrator: LangevinSplittingDynamicsMove")
print("Sampler: ReplicaExchangeSampler")

lsd_move = openmmtools.mcmc.LangevinSplittingDynamicsMove(timestep=timeStep, collision_rate=1.0/unit.picoseconds, n_steps=stepsPerIteration)
for k in range(nstates):
	sampler_state = sampler_states[k]
	thermodynamic_state = thermodynamic_states[k]
	integrator=mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
	context = thermodynamic_state.create_context(integrator)
	system = context.getSystem()
	for force in system.getForces():
		if isinstance(force, mm.CustomBondForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.HarmonicBondForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.HarmonicAngleForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.PeriodicTorsionForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomAngleForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.NonbondedForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomNonbondedForce):
			force.updateParametersInContext(context)
		elif isinstance(force, mm.CustomTorsionForce):
			force.updateParametersInContext(context)
	sampler_state.apply_to_context(context)
	initial_energy = thermodynamic_state.reduced_potential(context)
	print("Sampler state {}: initial energy {:8.3f}kT".format(k, initial_energy))
	mm.LocalEnergyMinimizer.minimize(context)
	sampler_state.update_from_context(context)
	final_energy = thermodynamic_state.reduced_potential(context)
	print("Sampler state {}: final energy {:8.3f}kT".format(k, final_energy))
	del context

if runflag == 'run':
	repex_simulation = ReplicaExchangeSampler(mcmc_moves=lsd_move, number_of_iterations=productionIterations)

tmp_dir = './trajectory/'
storage = os.path.join(tmp_dir, 'solvent.nc')
reporter = MultiStateReporter(storage, checkpoint_interval=iterationsPerCheckpoint)
if runflag != 'run':
    repex_simulation = ReplicaExchangeSampler.from_storage(reporter)
else:
    repex_simulation.create(thermodynamic_states, sampler_states, reporter)
    repex_simulation.equilibrate(equilibrationIterations)
if runflag == 'recover' or runflag == 'run':
        repex_simulation.run()
elif runflag == 'extend':
        repex_simulation.extend(extendIterations)

#will add all iterations even if coming from a previous restart
all_iters = repex_simulation.iteration
print('All iterations = {}'.format(all_iters))

analyzer = ReplicaExchangeAnalyzer(reporter)
iterations_to_analyze=int(all_iters/50)
print('Iterations to analyze = {}'.format(iterations_to_analyze))
for i in range(1, iterations_to_analyze+1):
	samples_to_analyze=i*50
	analyzer.max_n_iterations = samples_to_analyze
	Delta_f_ij, dDelta_f_ij = analyzer.get_free_energy()
	print("Relative free energy change in {0} = {1} +- {2}"
		.format('solvent', Delta_f_ij[0, nstates - 1]*kTtokcal, dDelta_f_ij[0, nstates - 1]*kTtokcal))

[matrix,eigenvalues,ineff]=analyzer.generate_mixing_statistics()
print("Mixing Stats")
print(matrix)
print(eigenvalues)
print(ineff)
