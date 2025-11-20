
function dbdPdP = op_D_DP_DP_gradv_Br_bot_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, dCdP, d2CdP2, mat, patch_list)
    
    nParams = size(dCdP, 3);

    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list

        [ind, vals] = op_D_DP_DP_gradv_Br_bot(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:), mat);

        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end

    dbdPdP = sptensor (cell2mat(indices), cell2mat(values), [spv.ndof, nParams, nParams]);
end
