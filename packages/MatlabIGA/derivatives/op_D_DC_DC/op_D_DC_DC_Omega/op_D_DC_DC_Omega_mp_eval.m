
function dVdCdC = op_D_DC_DC_Omega_mp_eval(spg, spgEval, msh, mshEval, patch_list)

    if (spg.npatch ~= msh.npatch)
        error ('op_D_DC_DC_Omega_mp: the number of patches does not coincide')
    end
    
    inds = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    
    for iPatch = patch_list
    
        [ind, vs] = op_D_DC_DC_Omega(spgEval{iPatch}, mshEval{iPatch});
    
        ind(:,1) = spg.gnum{iPatch}(ind(:,1));
        % ind(:,2) = ind(:,2);
        ind(:,3) = spg.gnum{iPatch}(ind(:,3));
        % ind(:,4) = ind(:,4);
    
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_DC_Omega_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iPatch};
        end
    
        inds{iPatch} = ind;
        vals{iPatch} = vs;
    end
    
    dVdCdC = sptensor (cell2mat(inds), cell2mat(vals), [spg.ndof, msh.ndim, spg.ndof, msh.ndim]);

end
