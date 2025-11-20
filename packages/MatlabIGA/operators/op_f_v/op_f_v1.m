% OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
%
%   rhs = op_f_v (spv, msh, func);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   func:  function handle
%

function varargout = op_f_v1 (spv, msh, func)
    if ~exist("func", "var")
        func = @(x,y,z) repmat(reshape(ones(size(x)), [1, size(x)]), spv.ncomp, 1, 1);
    end
    for idim = 1:msh.rdim
        x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
    end
    coeff = func(x{:});  % [ncomp, nqn, nel]
    coeff = reshape (coeff, spv.ncomp, msh.nqn, 1, msh.nel);
    
    shpv  = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel); % [ncomp, nqn, nv, nel]
    
    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,2]);
    
    values = reshape(sum(sum(shpv.*coeff, 1) .* jacdet, 2), [], 1);
    
    rows = reshape(spv.connectivity, [], 1);
    cols = ones(size(rows));
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, 1);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_f_v1: wrong number of output arguments')
    end
end

