classdef MAT_416steel < MAT_Magnetic & MAT_Thermal

    methods (Access = public)
        function obj = MAT_416steel()
            obj.Mur = 440;
            obj.Sigma = 1.03e7;
            obj.Rho = 7874;
            obj.IsLinearMAG = false;
            obj.B = [0.000000, 0.200400, 0.425600, 0.761000, 1.097000, 1.233000, 1.335000, 1.460000, 1.590000, 1.690000, 1.724232, 1.740000];
            obj.H = [0.000000, 318.310000, 477.465000, 795.775000, 1591.550000, 2387.325000, 3978.875000, 7957.750000, 15915.500000, 31831.000000, 44456.280000, 55704.250000];
            obj.fitHBspline();
        end
    end
end


