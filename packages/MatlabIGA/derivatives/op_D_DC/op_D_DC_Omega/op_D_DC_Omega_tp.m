function varargout = op_D_DC_Omega_tp(space1, msh)
    
    for idim = 1:msh.ndim
        size1 = size (space1.sp_univ(idim).connectivity);
        if (size1(2) ~= msh.nel_dir(idim))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    A = sparse(space1.ndof, msh.ndim);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);
    
        A = A + op_D_DC_Omega(sp_col, msh_col);
    end
    
    if (nargout == 1)
        varargout{1} = A;
    elseif (nargout == 3)
        [rows, dims, vals] = find (A);
        varargout{1} = rows;
        varargout{2} = dims;
        varargout{3} = vals;
    end
end