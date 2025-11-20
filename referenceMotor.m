% SPDX-License-Identifier: GPL-3.0

function [] = referenceMotor(simParams)

% motor1Az - Simulates and evaluates the magnetic potential of a parameterized
% permanent magnet synchronous motor (IPMSM) using isogeometric analysis.
%
% Syntax:
%   x = referenceMotor(simParams)
%
% INPUT:
%
%  simParams: a structure with data of the simulation run. It contains the fields:
%   iter       - Integer index for the simulation run (used in file naming).
%   mh      - Magnet height (geometry parameter).
%   mw      - Magnet width (geometry parameter).
%   mag     - Magnet distance to airgap (geometry parameter).
%   rot     - Rotor rotation angle in degrees.
%   airgap  - Boolean flag (true/false). If true, extract solution only in the airgap.
%   predictionsDir   - Directory of the prediction data.
%   predictionsFile    - Name of the prediction file.
%   errorFile   - Name where to write the errors in.
%
% Description:
%   This function sets up the geometry, materials, boundary conditions, and excitations
%   for a parameterized PMSM. It solves the magnetostatic problem using harmonic mortaring
%   and isogeometric analysis (IGA). It computes torque and field predictions errors.


addpath(genpath("packages"))

[solShrink, K, Rotor, Stator, Motor, Coupling, solFull] = solveMotor(simParams);

% Extract the fields from the data structures into local variables
data_names = fieldnames (simParams);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= simParams.(data_names{iopt});']);
end

if(airgap)
    AirGapDofs = [];
    AigGapPatches = 13:17;
    for iPatch = 1:numel(AigGapPatches)
        AirGapDofs = union(AirGapDofs, Rotor.Spaces.gnum{AigGapPatches(iPatch)});
    end
    x = Rotor.Solution.Values(AirGapDofs);
    Torque = Coupling.calcTorqueBrBtRt();
    
    % open prediction file
    pred_sol = readmatrix(predictionsDir + "/" + predictionsFile);
    K = load("MatrixK_gap.csv");
        
    Rotor.MagneticPotential = zeros(size(Rotor.Solution.Values));
    Rotor.MagneticPotential(AirGapDofs) = pred_sol;

    Torque_pred = Coupling.calcTorqueBrBtRt();
else
    Rind = Motor.getIndicesReduced(Rotor);
    Sind = Motor.getIndicesReduced(Stator);
    x = solShrink;
    Torque = Coupling.Solution.Torque;

    % open other prediction file
    pred_sol = readmatrix(predictionsDir + "/" + predictionsFile);
      
    newSolution = solFull(Motor.K);
    newSolution(Rind) = pred_sol(1:size(Rind,1));
    newSolution(Sind) = pred_sol(size(Rind,1)+1:end);
    newSolution = Motor.reconstructSolution(0, newSolution);
    Motor.postprocess(0, newSolution);

    Torque_pred = Coupling.Solution.Torque;

end

% Compute error in Hcurl seminorm
diff = pred_sol - x;
rel_err_curl =sqrt((diff'* K * diff) /(x'* K * x));
% Compute relative error of torque
rel_err_torque = abs(Torque-Torque_pred)/ Torque;

err_data = [iter, rel_err_curl, rel_err_torque];

% Save to error file
if isfile(errorFile)
    if iter==1
        msg = 'This error file already exists. Are you sure you want to write into an existing file?';
        error(msg)
    end
    fid = fopen(errorFile,'a');
else
    fid = fopen(errorFile, 'w');
    fprintf(fid, 'iter, curl_err, torque_err\n');
end

fprintf(fid, '%i,%f,%f\n', err_data);
fclose(fid);
end
