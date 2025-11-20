

function dVdPdP = op_D_DP_DP_Omega_mp_eval(spg, spgEval, mshEval, dCdP, d2CdP2, patch_list)
    nParams = size(dCdP, 3);
    dVdPdP = zeros(nParams, nParams);
    
    for iPatch = patch_list   

        values = op_D_DP_DP_Omega(spgEval{iPatch}, mshEval{iPatch}, dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:));

        dVdPdP = dVdPdP + values;
    end
end
