function [ControlPoints, all_params] = getC(geometryFunction, parameters, spaceGeo)
    parameters.draw_geometry = false;
    [srf, ~, all_params] = geometryFunction(parameters);

    for iPatch  = 1:spaceGeo.npatch
        ind_loc = spaceGeo.gnum{iPatch};
        ControlPoints(ind_loc, 1) = reshape(srf(iPatch).coefs(1, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 2) = reshape(srf(iPatch).coefs(2, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 3) = reshape(srf(iPatch).coefs(3, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
        ControlPoints(ind_loc, 4) = reshape(srf(iPatch).coefs(4, :, :, :), [], 1);
    end
end
