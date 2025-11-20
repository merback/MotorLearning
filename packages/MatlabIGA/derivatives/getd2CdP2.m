function d2CdP2 = getd2CdP2(geometryFunction, parameters, spaceGeo, Steps)
    nParameters = numel(fieldnames(parameters));
    if ~exist("Steps", "var")
        Steps = ones(nParameters, 1)*1e-6;
    end

    [dCdP0, all_params] = getdCdP(geometryFunction, parameters, spaceGeo, Steps);
    d2CdP2 = zeros([size(dCdP0), nParameters]);
    
    % Check, which parameters change the geometry first
    parameterNames = fieldnames(parameters);
    [~, indNonZero] = intersect(parameterNames, fieldnames(all_params), "stable");
    % indNonZero = 1:nParameters;

    for iParam = reshape(indNonZero, 1, [])
        Step = Steps(iParam);
        opts1 = parameters;
        opts1.(parameterNames{iParam}) = opts1.(parameterNames{iParam})+Step;
        dCdP1 = getdCdP(geometryFunction, opts1, spaceGeo, Steps);

        d2CdP2(:,:,:,iParam) = (dCdP1-dCdP0)/(Step);

        % opts2 = parameters;
        % opts2.(parameterNames{iParam}) = opts2.(parameterNames{iParam})-Step;
        % dCdP2 = getdCdP(geometryFunction, opts2, spaceGeo, Steps);
        % d2CdP2(:,:,:,iParam) = (dCdP1-dCdP2)/(2*Step);
    end
end