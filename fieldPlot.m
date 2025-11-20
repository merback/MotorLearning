% SPDX-License-Identifier: GPL-3.0

function [] = fieldPlot()

% fieldPlot - Reproduces the field plot in Fig. 5 from 10.1109/TMAG.2025.3629546
%
% Requirements:
%   - External data folders for predictions must exist.

addpath(genpath("packages"))
% Get the simulation parameters
filename = sprintf("samples_full_testing.csv");
data = readmatrix(filename);
i = 121;

simParams.iter = data(i,1);
simParams.mh = data(i, 2);
simParams.mw = data(i, 3);
simParams.magn = data(i, 4);
simParams.rot = data(i, 5);
simParams.airgap = false;
simParams.lw1 = data(i, 6);
simParams.lw2 = data(i, 7);
simParams.theta1 = data(i, 8);
simParams.theta2 = data(i, 9);
simParams.i1 = data(i,10);
simParams.i2 = data(i,11);
simParams.i3 = data(i, 12);
simParams.phi1 = data(i,13);
simParams.phi2 = data(i,14);
simParams.phi3 = data(i,15);

simParams.computeK = false;

[~, ~, Rotor, Stator, Motor, ~, solFull] = solveMotor(simParams);

Rind = Motor.getIndicesReduced(Rotor);
Sind = Motor.getIndicesReduced(Stator);

Motor.plotMagneticFluxDensity(0,2);
colormap(viridis)
colorbar
exportgraphics(gcf, "Fieldplot_groundtruth.pdf", "ContentType", "vector")

filename = sprintf('Predictions_full_testing/full_testing_%d.csv', i);
pred_sol = readmatrix(filename);
newSolution = solFull(Motor.K);
newSolution(Rind) = pred_sol(1:size(Rind,1));
newSolution(Sind) = pred_sol(size(Rind,1)+1:end);
newSolutionF = Motor.reconstructSolution(0, newSolution);
Motor.postprocess(0, newSolutionF);
cols = [0,2];
Motor.plotMagneticFluxDensity(0,2);
colormap(viridis)
clim(cols)
colorbar
exportgraphics(gcf, "Fieldplot_prediction.pdf", "ContentType", "vector");
        
% Plot absolute error
hold on
newSolution_err = Motor.reconstructSolution(0, solFull(Motor.K)-newSolution);
Motor.postprocess(0, newSolution_err);
Motor.plotMagneticFluxDensity(0,0.1);
colormap(viridis)
colorbar
exportgraphics(gcf, "Fieldplot_abserr.pdf", "ContentType", "vector")

