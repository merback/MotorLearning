
function dbdPdP = op_D_DP_DP_gradv_Br_bot_times_z_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, z, dCdP, d2CdP2, mat, patch_list)
    
    nParams = size(dCdP, 3);
    dbdPdP = zeros(nParams, nParams);
    
    for iPatch = patch_list

        vals = op_D_DP_DP_gradv_Br_bot_times_z(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, z(spv.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:), mat);

        dbdPdP = dbdPdP + vals;
    end

end
