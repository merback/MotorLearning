classdef MAT_Magnetic_eff < MAT_Magnetic

    methods (Access = public)
        function obj = MAT_Magnetic_eff()
            obj@MAT_Magnetic();
        end

        function mu = getMuNonlinear(obj, B)
            B(B<=1e-8) = 1e-8; % avoid zero division
            mu = B./ppval(obj.HBsplineEFF, B);
            % correct B values higher than measurement (assume behavior like vacuum)
            mu(B>obj.Beff(end)) = B(B>obj.Beff(end))./((B(B>obj.Beff(end))-obj.Beff(end))./obj.Mu0 + obj.H(end));
        end
    end    
end