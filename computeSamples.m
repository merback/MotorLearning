% SPDX-License-Identifier: GPL-3.0

function [] = computeSamples(IsGap)


    % Load the CSV file
    format long g

    if IsGap
        name = "gap";
    else
        name = "full";
    end
    % Loop through testing, training and validation data
    for type = ["testing", "validation", "training"]
        
        filename = sprintf("samples_%s_%s.csv", name, type);
        data = readmatrix(filename);
        
        simParams.airgap = IsGap;
        simParams.computeK = false;
        
        % Loop through each row
        for i =1:size(data, 1)            
            simParams.iter = data(i,1);
            simParams.mh = data(i, 2);
            simParams.mw = data(i, 3);
            simParams.magn = data(i, 4);
            simParams.rot = data(i, 5);
            simParams.airgap = IsGap;

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
         
            % Call function
            x = solveMotor(simParams);
               
            % Create folder and save the output to a csv file
            if simParams.iter == 1
                dir_name = sprintf("sol_%s_%s", name, type); 
                mkdir(dir_name); % Create directory if it doesn't exist
            end
            str = dir_name + sprintf("/sol_%s_%s_", name, type) + simParams.iter + ".csv";
            writematrix(x, str);
        end
    end

    disp('All samples processed.');
    disp('Compute the (mean) stiffness matrix now.');
    

    % Compute (mean) stiffness matrix
    clear simParams
    simParams.airgap = IsGap;
    simParams.computeK = true;
    simParams.mh = 7e-3;
    simParams.mw = 20e-3;
    simParams.magn = 10e-3;
    simParams.rot = 10;
    simParams.lw1 = 2e-3;
    simParams.lw2 = 3e-3;
    if IsGap
        simParams.theta1 = rad2deg(0.37);
        simParams.theta2 = rad2deg(0.00558);
    else
        simParams.theta1 = 20;
        simParams.theta2 = 0.35;
    end
    simParams.i1 = 10;
    simParams.i2 = 10;
    simParams.i3 = 10;
    simParams.phi1 = 0;
    simParams.phi2 = 0;
    simParams.phi3 = 0;

    [~, K,~,~,~,~,~] = solveMotor(simParams);
    % Save the computed stiffness matrix to a CSV file
    writematrix(full(K), sprintf('MatrixK_%s.csv', name));
    disp('Stiffness matrix computed.');
end
