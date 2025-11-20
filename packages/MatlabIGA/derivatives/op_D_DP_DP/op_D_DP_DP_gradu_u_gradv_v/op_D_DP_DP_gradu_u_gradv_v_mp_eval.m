
function dKdPdP = op_D_DP_DP_gradu_u_gradv_v_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, v, dCdP, d2CdP2, patch_list)
    
    nParams = size(dCdP, 3);
    dKdPdP = zeros(nParams, nParams);
    
    for iPatch = patch_list

        vals = op_D_DP_DP_gradu_u_gradv_v(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), v(spv.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:));

        dKdPdP = dKdPdP + vals;
    end
end
