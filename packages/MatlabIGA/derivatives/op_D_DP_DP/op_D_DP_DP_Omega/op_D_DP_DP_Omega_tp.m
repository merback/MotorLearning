
function varargout = op_D_DC_DC_Omega_tp(spaceGeo, msh)

    for idim = 1:msh.ndim
        sizeGeo = size (spaceGeo.sp_univ(idim).connectivity);
        if sizeGeo(2) ~= msh.nel_dir(idim)
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    dVdCdC = sptensor([],[],[spaceGeo.ndof, msh.ndim, spaceGeo.ndof, msh.ndim]);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        % spGeo_col = sp_evaluate_col_param(spaceGeo, msh_col, 'value', true, 'gradient', true);
        spGeo_col = sp_evaluate_col(spaceGeo, msh_col, 'value', true, 'gradient', true);
    
        dVdCdC = dVdCdC + op_D_DC_DC_Omega(spGeo_col, msh_col);
    end
    
    if (nargout == 1)
        varargout{1} = dVdCdC;
    elseif (nargout == 2)
        [indices, vals] = find (dVdCdC);
        varargout{1} = indices;
        varargout{2} = vals;
    end
end