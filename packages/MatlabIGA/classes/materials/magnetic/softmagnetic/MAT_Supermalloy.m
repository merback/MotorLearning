classdef MAT_Supermalloy < MAT_Magnetic & MAT_Thermal

    methods (Access = public)
        function obj = MAT_Supermalloy()
            obj.Mur = 4.5e5;
            obj.Sigma = 1.03e7;
            obj.Rho = 7874;
            obj.IsLinearMAG = false;
            obj.B = [0.000000, 0.450001, 0.575002, 0.625003, 0.650004, 0.680006, 0.700008, 0.725673, 0.744506, 0.760030, 0.780460, 0.785100];
            obj.H = [0.000000, 0.795775, 1.591550, 2.387325, 3.183100, 4.774649, 6.366200, 10.138480, 15.685950, 23.873250, 48.000000, 79.577500];
            obj.fitHBspline();
        end
    end
end