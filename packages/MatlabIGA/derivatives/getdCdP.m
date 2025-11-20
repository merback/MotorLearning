function [dCdP, all_params] = getdCdP(geometryFunction, parameters, spaceGeo, Steps)
    nParameters = numel(fieldnames(parameters));
    if ~exist("Steps", "var")
        Steps = ones(nParameters, 1)*1e-6;
    end

    [C0, all_params] = getC(geometryFunction, parameters, spaceGeo);
    dCdP = zeros([size(C0), nParameters]);

    % Check, which parameters change the geometry first
    parameterNames = fieldnames(parameters);
    [~, indNonZero] = intersect(parameterNames, fieldnames(all_params), "stable");
    % indNonZero = 1:nParameters;
    
    for iParam = reshape(indNonZero, 1, [])
        Step = Steps(iParam);
        opts1 = parameters;
        opts1.draw_geometry = false;
        opts1.(parameterNames{iParam}) = opts1.(parameterNames{iParam})+Step;
        C1 = getC(geometryFunction, opts1, spaceGeo);

        dCdP(:, :, iParam) = (C1 - C0)/(Step);

        % opts2 = parameters;
        % opts2.draw_geometry = false;
        % opts2.(parameterNames{iParam}) = opts2.(parameterNames{iParam})-Step;
        % C2 = getC(geometryFunction, opts2, spaceGeo);
        % 
        % dCdP(:, :, iParam) = (C1 - C2)/(2*Step);
    end
end
