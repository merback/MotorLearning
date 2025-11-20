% SPDX-License-Identifier: GPL-3.0

function [] = computeErrors(IsGap, paper)


    % Load the CSV file
    format long g

    if IsGap
        name = "gap";
    else
        name = "full";
    end

    % Loop through testing, training and validation data
    for type = ["validation", "testing","training"]
        
        filename = sprintf("samples_%s_%s.csv", name, type);
        data = readmatrix(filename);
        
        simParams.airgap = IsGap;
        simParams.predictionsDir = sprintf("Predictions_%s_%s/",name, type);
        simParams.errorFile = sprintf("Error_%s_%s.csv", name, type);
        if paper
            simParams.predictionsDir = sprintf("Predictions_paper_%s_%s/", name, type);
            simParams.errorFile = sprintf("Error_paper_%s_%s.csv", name, type);
        end
        
        % Loop through each row
        for i =1:size(data, 1)            
            simParams.iter = data(i,1);
            simParams.mh = data(i, 2);
            simParams.mw = data(i, 3);
            simParams.magn = data(i, 4);
            simParams.rot = data(i, 5);
            simParams.airgap = IsGap;
            simParams.predictionsFile = sprintf("%s_%s_%i", name, type, i);

            if ~IsGap
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
            else
                simParams.lw1 = 2e-3;
                simParams.lw2 = 3e-3;
                simParams.theta1 = rad2deg(0.37);
                simParams.theta2 = rad2deg(0.00558);
                simParams.i1 = 10;
                simParams.i2 = 10;
                simParams.i3 = 10;
                simParams.phi1 = 0;
                simParams.phi2 = 0;
                simParams.phi3 = 0;
            end
            
            simParams.computeK = true;
            % Call function
            referenceMotor(simParams);
            
        end
        fprintf("Errors successfully compute for %s data. ",type)
    end

    disp('All samples processed. Errors computed.');
end
