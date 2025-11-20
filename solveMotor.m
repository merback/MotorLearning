% SPDX-License-Identifier: GPL-3.0

function [x, K, Rotor, Stator, Motor, Coupling, sol] = solveMotor(simParams)

% motor1Az - Simulates and evaluates the magnetic potential of a parameterized
% permanent magnet synchronous motor (IPMSM) using isogeometric analysis.
%
% Syntax:
%   [x, K] = solveMotor(simParams)
%
% INPUT:
%
%  simParams: a structure with data of the simulation run. It contains the fields:
%   iter     - Integer index for the simulation run (used in file naming).
%   mh       - Magnet height (geometry parameter).
%   mw       - Magnet width (geometry parameter).
%   mag      - Magnet distance to airgap (geometry parameter).
%   rot      - Rotor rotation angle in degrees.
%   airgap   - Boolean flag (true/false). If true, extract solution only in the airgap.
%   computeK - Boolean flag (true/false). If true, compute the stiffness matrix with unit material
%   lw1
%   lw2
%   theta1
%   theta2
%   i1
%   i2
%   i3
%   phi1
%   phi2
%   phi3
%
% Outputs:
%   x       - Solution vector. Either the magnetic potential in the airgap or the full solution.
%   K       - Stiffness matrix.
%
% Description:
%   This function sets up the geometry, materials, boundary conditions, and excitations
%   for a parameterized PMSM. It solves the magnetostatic problem using harmonic mortaring
%   and isogeometric analysis (IGA) and returns the solution vector.


addpath(genpath("packages"))

% Extract the fields from the data structures into local variables
data_names = fieldnames (simParams);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= simParams.(data_names{iopt});']);
end

%% Rotor
% Geometry
rotationAngle = rot;
params.draw_geometry=false;
params.MW = mw;
params.MH = mh;
params.MAG = magn;
if ~airgap
    params.LW1 = lw1;
    params.LW2 = lw2;
    params.Theta1 = deg2rad(theta1);
    params.Theta2 = deg2rad(theta1 + theta2);
end

[srfRotor, patchesRotor, paramsRt] = GEO_IPMrotor_I_1(params);
[srfStator, patchesStator, paramsSt] = GEO_DISstator_1(params);

%%
Rotor = IGA_RegionMagnetic();
Rotor.importSurface(srfRotor);
Rotor.setProperties([5,5], [2,2], 0.1);
Rotor.calculatePlotPoints([5, 5]);
Stator = IGA_RegionMagnetic();
Stator.importSurface(srfStator);
Stator.setProperties([4,4], [2,2], 0.1);
Stator.calculatePlotPoints([5,5]);

%% Materials
Air = MAT_Air();
Iron = MAT_M27();
Copper = MAT_Copper();
Magnet = MAT_NdFeB_Br14(pi/2);

Rotor.resetMaterials();
Rotor.setMaterial(Air, patchesRotor.Air);
Rotor.setMaterial(Iron, patchesRotor.Iron);
Rotor.setMaterial(Magnet, patchesRotor.Magnets);

Stator.resetMaterials();
Stator.setMaterial(Air, patchesStator.Air);
Stator.setMaterial(Iron, patchesStator.Iron);
Stator.setMaterial(Copper, patchesStator.Windings);

%%
Rotor.resetBoundaryConditions()
Rotor.addBoundaryCondition(IGA_BCdirichlet(Rotor, 4));
Rotor.addBoundaryCondition(IGA_BCantiPeriodic(Rotor, [5,2,11], [3,1,6]));

%%
Stator.resetBoundaryConditions()
Stator.addBoundaryCondition(IGA_BCdirichlet(Stator, [8,12,16,20,24,32]));
Stator.addBoundaryCondition(IGA_BCantiPeriodic(Stator, [27,29,30,31], [1,5,6,7]));

%%
phaseDefU(1).Patches = [8, 10, 13, 14];        phaseDefU(1).Slot = 1;
phaseDefU(2).Patches = [7, 9, 11, 12];         phaseDefU(2).Slot = 1;
phaseDefU(3).Patches = [26, 28, 30, 31];       phaseDefU(3).Slot = 1;
phaseDefU(4).Patches = [103, 105, 108, 109];   phaseDefU(4).Slot = -1;

phaseDefW(1).Patches = [27, 29, 32, 33];       phaseDefW(1).Slot = -1;
phaseDefW(2).Patches = [46, 48, 51, 52];       phaseDefW(2).Slot = -1;
phaseDefW(3).Patches = [45, 47, 49, 50];       phaseDefW(3).Slot = -1;
phaseDefW(4).Patches = [64, 66, 68, 69];       phaseDefW(4).Slot = -1;

phaseDefV(1).Patches = [65,67,70,71];          phaseDefV(1).Slot = 1;
phaseDefV(2).Patches = [83,85,87,88];          phaseDefV(2).Slot = 1;
phaseDefV(3).Patches = [84,86,89,90];          phaseDefV(3).Slot = 1;
phaseDefV(4).Patches = [102,104,106,107];      phaseDefV(4).Slot = 1;

NumberWindings = 12;


pp = 3;
Iu = i1; 
Iv = i2; 
Iw = i3; 
phase0 = 0;
PhaseU = deg2rad(phi1);
PhaseV = deg2rad(-120+phi2);
PhaseW = deg2rad(120+phi3);

% synchronous speed: omega*t = theta*pp
PhaseU = IGA_EXCcurrentSine(Stator, phaseDefU, Iu, pp, PhaseU, phase0, NumberWindings);
PhaseU.Color = TUDa_getColor("9c");
PhaseV = IGA_EXCcurrentSine(Stator, phaseDefV, Iv, pp, PhaseV, phase0, NumberWindings);
PhaseV.Color = TUDa_getColor("4c");
PhaseW = IGA_EXCcurrentSine(Stator, phaseDefW, Iw, pp, PhaseW, phase0, NumberWindings);
PhaseW.Color = TUDa_getColor("1c");


Stator.resetExcitations();
Stator.addExcitation(PhaseU);
Stator.addExcitation(PhaseV);
Stator.addExcitation(PhaseW);

%% Coupling
Coupling = IGA_HarmonicMortaring(Rotor, [7,8,9,10,12], Stator, [2,3,4,9,10,11,13,14,15,17,18,19,21,22,23,25,26,28]);
nharm = 25; 
mult = 6;
CosValues = 3:mult:mult*nharm/2;
SinValues = 3:mult:mult*nharm/2;
Coupling.setHarmonics(SinValues, CosValues);
%% Coupled System
Motor = IGA_MotorClass(Rotor, Stator);
Motor.setCoupling(Coupling);

%% Solve
currentAngle = deg2rad(rotationAngle);
Motor.setRotationAngle(currentAngle, "rad");
sol = Motor.solveStatic(currentAngle);

if(airgap)
    AirGapDofs = [];
    AigGapPatches = 13:17;
    for iPatch = 1:numel(AigGapPatches)
        AirGapDofs = union(AirGapDofs, Rotor.Spaces.gnum{AigGapPatches(iPatch)});
    end
    x = Rotor.Solution.Values(AirGapDofs);
else
    solShrink = sol(Motor.K);
    Rind = Motor.getIndicesReduced(Rotor);
    Sind = Motor.getIndicesReduced(Stator);
    x = [solShrink(Rind); solShrink(Sind)];
end

if(computeK)
    K_r = op_gradu_gradv_mp(Rotor.Spaces, Rotor.Spaces, Rotor.Meshes,@(x,y) ones(size(x)));
    if(airgap)
        K = K_r(AirGapDofs, AirGapDofs);
    else
        K_s = op_gradu_gradv_mp(Stator.Spaces, Stator.Spaces, Stator.Meshes, @(x,y) ones(size(x)));
        
        %Assemble corresponding K Matrix
        K = zeros(size(x,1),size(x,1));
        K(1:size(Rind,1), 1:size(Rind,1)) = K_r(Rotor.K,Rotor.K);
        K((size(Rind,1)+1):end, (size(Rind,1)+1):end) = K_s(Stator.K, Stator.K);
    end
else
    K = 0;
end
