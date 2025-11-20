classdef MAT_SoftIron < MAT_Magnetic & MAT_Thermal

    methods (Access = public)
        function obj = MAT_SoftIron()
            obj.Mur = 1200;
            obj.Sigma = 1.03e7;
            obj.Rho = 7874;
            obj.IsLinearMAG = false;
            obj.H = [0, 663.146/2, 663.146, 1067.5,1705.23, 2463.11, 3841.67, 5425.74, 7957.75, 12298.3, 20462.8, 32169.6, 61213.4, 111408, 188487.757, 267930.364, 347507.836];
            obj.B = [0, 0.5, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4];
            obj.fitHBspline();
        end
    end
end