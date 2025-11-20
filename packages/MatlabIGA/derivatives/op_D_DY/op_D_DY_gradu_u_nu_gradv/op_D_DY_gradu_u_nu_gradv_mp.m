function A = op_D_DY_gradu_u_nu_gradv_mp(spu, spv, msh, u, material, patch_list)
    
    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);

    for iptc = patch_list
        [r, c, vals] = op_D_DY_gradu_u_nu_gradv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, u(spu.gnum{iptc}), material);
        if isempty(r)
            continue
        end

        rows{iptc} = spv.gnum{iptc}(r);
        cols{iptc} = spu.gnum{iptc}(c);
        
        if (~isempty (spu.dofs_ornt))
            vals = spu.dofs_ornt{iptc}(rs)' .* vals;
        end
        if (~isempty (spv.dofs_ornt))
            vals = vals .* spv.dofs_ornt{iptc}(cs)';
        end
    
        values{iptc} = vals;
    end
    
    A = sparse(cell2mat(rows), cell2mat(cols), cell2mat(values), spu.ndof, spv.ndof);
end
