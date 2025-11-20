classdef IGA_Region < DAE_Element

    properties
        % quality of computation
        FormFunctionDegree ;                % degree of form functions 1-linear, 2-quadratic, 3-cubic, ...
        DegreeQuadrature;                   % degree of quadrature, usually use FormFunctionDegree
        SubdivisionsPatches;                % subdivision for each patch

        % domain definitions
        Length;                             % axial length
        Surface;                            % NURBS representation of single patches
        Geometry;                           % parametrization of the mapping to the physical domain and derivatives
        Boundaries;                         % information about the patches with an outer boundary
        Interfaces;                         % information about the connections between patches
        BoundaryInterfaces;                 % boundary information
        ControlPoints;                      % control points and weights

        % computation of spaces and meshes
        Spaces, SpacesEval;                 % space stores the shape functions of discrete spaces
        SpacesGeo, SpacesGeoEval;           % space of mapping functions for the geometry
        Meshes, MeshesEval;                 % mesh stores partition of domain and quadrature rule
        NumberPatches;                      % patch counter

        Materials;

        EvaluationPoints, EvaluationPatches, EvaluationShapeValues;
        EvaluationValues;

        % plotting
        FigureGeometry, FigureGeometryPatches, FigureGeometryBoundaries, FigureGeometryInterfaces;
        PlotLimitsX = [-Inf, Inf];
        PlotLimitsY = [-Inf, Inf];
        PlotResolution = [1, 1];
        RotationAngle = 0;      % current rotation angle of geometry
        PlotAngle = 0;          % angle when geometry was last plotted, if ~=RotationAngle, update plot points
    end

    methods (Abstract)
        generateSpaces(obj)
    end

    methods
        function obj = IGA_Region()
            obj.resetMaterials();
            obj.resetBoundaryConditions();
            obj.resetExcitations();
            obj.FigureGeometry = gobjects(1);
        end

        function setProperties(obj, subdivisions, degrees, length, name)
            obj.SubdivisionsPatches = subdivisions;
            obj.FormFunctionDegree = degrees;
            obj.DegreeQuadrature = degrees+1;
            obj.Length = length;
            if ~exist('name', 'var')
                name = inputname(1);
            end
            obj.Name = name;
            obj.generateMeshes();
            obj.generateSpaces();
            obj.initializeFigures();
        end

        function importSurface(obj, surface)
            % warning("The warning nrbderiv:SecondDerivative is turned off");
            warning('off', 'nrbderiv:SecondDerivative')
            obj.Surface = surface;
            [obj.Geometry, obj.Boundaries, obj.Interfaces, ~, obj.BoundaryInterfaces] = mp_geo_load(obj.Surface, 1e-8);
            obj.NumberPatches = numel(obj.Geometry);
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).PlotColor = [1,1,1];
                obj.Geometry(iPatch).PlotAlpha = 1;
                obj.Geometry(iPatch).Material = [];
            end
        end

        function generateMeshes(obj)
            progressbar([obj.Name ': Generating meshes: '])
            % pre-evaluated spaces and meshes
            obj.MeshesEval = cell (1, obj.NumberPatches);
            obj.SpacesGeoEval = cell (1, obj.NumberPatches);
            % temporary cell arrays
            meshes = cell (1, obj.NumberPatches);
            spacesGeo = cell (1, obj.NumberPatches);

            % create geometry information
            for iPatch = 1:obj.NumberPatches
                % physical mapping and quadrature
                [~, zeta] = kntrefine (obj.Geometry(iPatch).nurbs.knots, obj.SubdivisionsPatches-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
                rule      = msh_gauss_nodes (obj.DegreeQuadrature);
                [qn, qw]  = msh_set_quad_nodes (zeta, rule);
                meshes{iPatch} = msh_cartesian (zeta, qn, qw, obj.Geometry(iPatch));
                spacesGeo{iPatch} = sp_nurbs(obj.Geometry(iPatch).nurbs, meshes{iPatch});

                % mesh and mapping evaluation
                obj.MeshesEval{iPatch} = msh_precompute(meshes{iPatch});
                obj.SpacesGeoEval{iPatch}  = sp_precompute(spacesGeo{iPatch}, meshes{iPatch}, 'value', true, 'gradient', true);

                progressbar(iPatch/obj.NumberPatches);
            end

            obj.Meshes = msh_multipatch(meshes, obj.Boundaries);
            obj.SpacesGeo = sp_multipatch(spacesGeo, obj.Meshes, obj.Interfaces, obj.BoundaryInterfaces);

            %  update control points
            obj.ControlPoints = zeros(obj.SpacesGeo.ndof, 4);
            globalNumbering = obj.SpacesGeo.gnum;   % fetch once to save time. Matlab does some strange copying otherwise
            for iPatch = 1:obj.NumberPatches
                localIndex = globalNumbering{iPatch};
                % localIndex = obj.SpacesGeo.gnum{iPatch};  % instead of this, it takes long time for some reason
                obj.ControlPoints(localIndex, 1) = reshape(obj.Geometry(iPatch).nurbs.coefs(1, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 2) = reshape(obj.Geometry(iPatch).nurbs.coefs(2, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 3) = reshape(obj.Geometry(iPatch).nurbs.coefs(3, :, :, :)./obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
                obj.ControlPoints(localIndex, 4) = reshape(obj.Geometry(iPatch).nurbs.coefs(4, :, :, :), [], 1);
            end

            obj.calculateAreaAndInertia();
            obj.calculatePlotPoints();

            progressbar('finished');
        end

        function resetMaterials(obj)
            obj.Materials = struct('Material', {}, 'Patches', {}, 'Mass', {}, 'Inertia', {});
        end

        function setMaterial(obj, material, patches)
            pos = numel(obj.Materials)+1;
            obj.Materials(pos).Material = material;
            obj.Materials(pos).Patches = patches;
            obj.Materials(pos).Mass = material.getRho()*sum(vertcat(obj.Geometry(patches).Area))*obj.Length;
            obj.Materials(pos).Inertia = material.getRho()*sum(vertcat(obj.Geometry(patches).Inertia))*obj.Length;
            for iPatch = patches
                obj.Geometry(iPatch).Material = material;
                obj.Geometry(iPatch).PlotColor = material.getPlotColor();
                obj.Geometry(iPatch).PlotAlpha = material.getPlotAlpha();
            end
            % Check if material was set for different patches already
            for iMat = 1:numel(obj.Materials)-1
                doublePatches = intersect(obj.Materials(iMat).Patches, patches);
                if ~isempty(doublePatches)
                    stayPatches = setdiff(obj.Materials(iMat).Patches, doublePatches);
                    obj.Materials(iMat).Patches = stayPatches;
                    obj.Materials(iMat).Mass = material.getRho()*sum(vertcat(obj.Geometry(stayPatches).Area))*obj.Length;
                    obj.Materials(iMat).Inertia = material.getRho()*sum(vertcat(obj.Geometry(stayPatches).Inertia))*obj.Length;
                    warning(['Materials for already specified patches ' num2str(doublePatches) ' were overwritten!']);
                end
            end
        end

        function dofs = getBoundaryDOFs(obj, boundaryNumbers)
            dofs = [];
            for bnd = boundaryNumbers
                localDOFS = obj.Spaces.boundary.gnum{bnd};
                globalDOFS = obj.Spaces.boundary.dofs(localDOFS);
                dofs = union(dofs, globalDOFS);
            end
        end

        % Sort DOFs such that the boundary discretization is sorted with the geometric sorting
        function dofs = getBoundaryDOFsSorted(obj, boundaryNumbers)
            dofs = [];
            for iBnd = 1:numel(boundaryNumbers)-1
                % first line element
                bndNr1 = boundaryNumbers(iBnd);
                patch1 = obj.Boundaries(bndNr1).patches;
                side1 = obj.Boundaries(bndNr1).faces;
                nrb1 = nrbextract(obj.Geometry(patch1).nurbs, side1);
                localDOFs1 = obj.Spaces.boundary.gnum{bndNr1};
                globalDOFs1 = obj.Spaces.boundary.dofs(localDOFs1);
                % second line element
                bndNr2 = boundaryNumbers(iBnd+1);
                patch2 = obj.Boundaries(bndNr2).patches;
                side2 = obj.Boundaries(bndNr2).faces;
                nrb2 = nrbextract(obj.Geometry(patch2).nurbs, side2);
                localDOFs2 = obj.Spaces.boundary.gnum{bndNr2};
                globalDOFs2 = obj.Spaces.boundary.dofs(localDOFs2);
                % check if mapping is reversed
                tol = 1e-10;
                if abs(vecnorm(nrbeval(nrb1, 1) - nrbeval(nrb2, 0))) < tol || abs(vecnorm(nrbeval(nrb1, 1) - nrbeval(nrb2, 1))) < tol
                    dofs = union(dofs, globalDOFs1, 'stable');
                elseif abs(vecnorm(nrbeval(nrb1, 0) - nrbeval(nrb2, 0))) < tol || abs(vecnorm(nrbeval(nrb1, 0) - nrbeval(nrb2, 1))) < tol
                    dofs = union(dofs, flip(globalDOFs1), 'stable');
                else
                    error("Boundaries are not matching or numbers not given sorted!")
                end
            end
            localDOFs1 = obj.Spaces.boundary.gnum{boundaryNumbers(end)};
            globalDOFs1 = obj.Spaces.boundary.dofs(localDOFs1);
            % Check last line segment
            if abs(vecnorm(nrbeval(nrb2, 0) - nrbeval(nrb1, 0))) < tol || abs(vecnorm(nrbeval(nrb2, 0) - nrbeval(nrb1, 1))) < tol
                dofs = union(dofs, globalDOFs1, 'stable');
            elseif abs(vecnorm(nrbeval(nrb2, 1) - nrbeval(nrb1, 0))) < tol || abs(vecnorm(nrbeval(nrb2, 1) - nrbeval(nrb1, 1))) < tol
                dofs = union(dofs, flip(globalDOFs1), 'stable');
            else
                error("Boundaries are not matching or numbers not given sorted!")
            end
        end

        function [leftDofs, rightDofs] = getBoundaryDOFsperiodic(obj, boundaryNumbersLeft, boundaryNumbersRight)
            assert(numel(boundaryNumbersLeft) == numel(boundaryNumbersRight), "Number of left and right sides must be the same!");
            leftDofs = [];
            rightDofs = [];
            for iBnd = 1:numel(boundaryNumbersRight)
                % extract left boundary
                bndNrLeft = boundaryNumbersLeft(iBnd);
                patchLeft = obj.Boundaries(bndNrLeft).patches;
                sideLeft = obj.Boundaries(bndNrLeft).faces;
                localDOFsLeft = obj.Spaces.boundary.gnum{bndNrLeft};
                globalDOFsLeft = obj.Spaces.boundary.dofs(localDOFsLeft);
                nrbLeft = nrbextract(obj.Geometry(patchLeft).nurbs, sideLeft);
                % extract right boundary
                bndNrRight = boundaryNumbersRight(iBnd);
                patchRight = obj.Boundaries(bndNrRight).patches;
                sideRight = obj.Boundaries(bndNrRight).faces;
                localDOFsRight = obj.Spaces.boundary.gnum{bndNrRight};
                globalDOFsRight = obj.Spaces.boundary.dofs(localDOFsRight);
                nrbRight = nrbextract(obj.Geometry(patchRight).nurbs, sideRight);
                % check if mapping is reversed
                tol = 1e-10;
                if abs(vecnorm(nrbeval(nrbRight, 0)) - vecnorm(nrbeval(nrbLeft, 0))) < tol
                    leftDofs = union(leftDofs, globalDOFsLeft, 'stable');
                    rightDofs = union(rightDofs, globalDOFsRight, 'stable');
                elseif abs(vecnorm(nrbeval(nrbRight, 0)) - vecnorm(nrbeval(nrbLeft, 1))) < tol
                    leftDofs = union(leftDofs, globalDOFsLeft, 'stable');
                    rightDofs = union(rightDofs, flip(globalDOFsRight), 'stable');
                else
                    error("Boundaries are not matching!")
                end
            end
        end

        %%%%%%%%%%%%%%%% misc
        function calculateAreaAndInertia(obj)
            % patches area and inertia
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).Area = sum(obj.MeshesEval{iPatch}.jacdet .* obj.MeshesEval{iPatch}.quad_weights, 'all');
                obj.Geometry(iPatch).Inertia = sum(obj.MeshesEval{iPatch}.jacdet .* obj.MeshesEval{iPatch}.quad_weights .* reshape((obj.MeshesEval{iPatch}.geo_map(1, :, :).^2 + obj.MeshesEval{iPatch}.geo_map(2, :, :).^2) , obj.MeshesEval{iPatch}.nqn, obj.MeshesEval{iPatch}.nel), 'all');
            end
            % interfaces lengths
            for iInt = 1:numel(obj.Interfaces)
                patch1 = obj.Interfaces(iInt).patch1;
                side1 = obj.Interfaces(iInt).side1;
                patch2 = obj.Interfaces(iInt).patch2;
                side2 = obj.Interfaces(iInt).side2;
                nrb = obj.Geometry(patch1).boundary(side1).nurbs;
                rule      = msh_gauss_nodes (nrb.order);
                [qn, qw]  = msh_set_quad_nodes (nrb.knots, rule);
                mesh = msh_cartesian(nrb.knots, qn, qw, geo_load(nrb));
                obj.Interfaces(iInt).length = op_Omega(msh_precompute(mesh));
                obj.Geometry(patch1).Lengths{side1} = obj.Interfaces(iInt).length;
                obj.Geometry(patch2).Lengths{side2} = obj.Interfaces(iInt).length;
            end
            % boundary lengths
            for iBnd = 1:numel(obj.Boundaries)
                patch = obj.Boundaries(iBnd).patches;
                side = obj.Boundaries(iBnd).faces;
                nrb = obj.Geometry(patch).boundary(side).nurbs;
                rule      = msh_gauss_nodes (nrb.order);
                [qn, qw]  = msh_set_quad_nodes (nrb.knots, rule);
                mesh = msh_cartesian(nrb.knots, qn, qw, geo_load(nrb));
                obj.Boundaries(iBnd).length = op_Omega(msh_precompute(mesh));
                obj.Geometry(patch).Lengths{side} = obj.Boundaries(iBnd).length;
            end
        end

        function calculatePlotPoints(obj, nsubs, patches)
            if exist("nsubs", "var")
                obj.PlotResolution = nsubs;
            end
            if ~exist("patches", "var")
                patches = 1:obj.NumberPatches;
            end

            for iPatch = patches
                pointsX = unique(kntrefine(obj.Geometry(iPatch).nurbs.knots{1}, obj.PlotResolution(1), obj.Geometry(iPatch).nurbs.order(1), obj.Geometry(iPatch).nurbs.order(1)-1));
                pointsY = unique(kntrefine(obj.Geometry(iPatch).nurbs.knots{2}, obj.PlotResolution(2), obj.Geometry(iPatch).nurbs.order(2), obj.Geometry(iPatch).nurbs.order(2)-1));
                plotPointsParam = {pointsX, pointsY};
                obj.Geometry(iPatch).PlotPointsParam = plotPointsParam;
                pnts = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, plotPointsParam, obj.PlotAngle);
                obj.Geometry(iPatch).PlotPoints = pnts(:, 1);
                % update interfaces
                for iInt = 1:numel(obj.Interfaces)
                    if obj.Interfaces(iInt).patch1 == iPatch
                        obj.Interfaces(iInt).PlotPoints = pnts(:, obj.Interfaces(iInt).side1+1);
                        obj.Interfaces(iInt).PlotColor = [0.5, 0.5, 0.5];
                    end
                end
                % update boundaries
                for iBnd = 1:numel(obj.Boundaries)
                    if obj.Boundaries(iBnd).patches == iPatch
                        obj.Boundaries(iBnd).PlotPoints = pnts(:, obj.Boundaries(iBnd).faces+1);
                        obj.Boundaries(iBnd).PlotColor = [0, 0, 0];
                    end
                end
            end
        end

        function updatePlotPoints(obj, rot)
            for iPatch = 1:obj.NumberPatches
                pnts = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsParam, rot);
                obj.Geometry(iPatch).PlotPoints = pnts(:, 1);
                % update interfaces
                for iInt = 1:numel(obj.Interfaces)
                    if obj.Interfaces(iInt).patch1 == iPatch
                        obj.Interfaces(iInt).PlotPoints = pnts(:, obj.Interfaces(iInt).side1+1);
                    end
                end
                % update boundaries
                for iBnd = 1:numel(obj.Boundaries)
                    if obj.Boundaries(iBnd).patches == iPatch
                        obj.Boundaries(iBnd).PlotPoints = pnts(:, obj.Boundaries(iBnd).faces+1);
                    end
                end
            end
        end

        function setEvaluationPoints(obj, evalPoints, evalPatches)
            obj.EvaluationPoints = evalPoints;
            obj.EvaluationPatches = evalPatches;
            obj.EvaluationShapeValues = shp_eval_mp(obj.Spaces, obj.Geometry, evalPoints, evalPatches, {'value', 'gradient'});
        end

        %%%%%%%%%%%%%%%% plots
        function plotGeometry(obj)
            if ~ishandle(obj.FigureGeometry)
                obj.FigureGeometry = figure('Name', 'Geometry Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
            else
                % figure(obj.FigureGeometry)        % focus current figure
            end
            hold on
            obj.plotPatches();
            obj.plotBoundaries(true);
            % set(obj.FigureGeometry.CurrentAxes, "XLim", obj.PlotLimitsX)
            % set(obj.FigureGeometry.CurrentAxes, "YLim", obj.PlotLimitsY)
            try
                axis equal
            catch
                warning('Something went wrong with axis equal');
            end
        end

        function plotPatches(obj)
            obj.checkRotationUpdate();
            for iMat = 1:numel(obj.Materials)
                for patch = obj.Materials(iMat).Patches
                    col = obj.Geometry(patch).PlotColor;
                    alph = obj.Geometry(patch).PlotAlpha;
                    if ishandle(obj.FigureGeometryPatches(patch))
                        set(obj.FigureGeometryPatches(patch), 'XData', obj.Geometry(patch).PlotPoints{1}, 'YData', obj.Geometry(patch).PlotPoints{2}, 'FaceColor', col, 'FaceAlpha', alph);
                    else
                        obj.FigureGeometryPatches(patch) = surf(obj.Geometry(patch).PlotPoints{1}, obj.Geometry(patch).PlotPoints{2}, zeros(size(obj.Geometry(patch).PlotPoints{1})), ...
                            "EdgeColor", "none", "FaceColor", col, "FaceAlpha", alph);
                    end
                end
            end
            % Check if all materials are defined
            allPatches = horzcat(obj.Materials.Patches);
            nonDefined = setdiff(1:obj.NumberPatches, allPatches);
            if ~isempty(nonDefined)
                disp("Patches " + num2str(nonDefined) + " were not set a material!")
            end
        end

        function plotKnotLines(obj)
            for iPatch = 1:obj.NumberPatches
                msh = obj.Meshes.msh_patch{iPatch};
                for knt = unique(obj.Spaces.sp_patch{iPatch}.knots{2})
                    knot_vals =  linspace(0,1,20);
                    xy = msh.map({knot_vals, knt});
                    line(xy(1,:), xy(2,:), 'Color', [1,1,1]*0.5);
                end
                % Y lines
                for knt = unique(obj.Spaces.sp_patch{iPatch}.knots{1})
                    knot_vals =  linspace(0,1,20);
                    xy = msh.map({knt, knot_vals});
                    line(xy(1,:), xy(2,:), 'Color', [1,1,1]*0.5);
                end
            end
        end

        function plotBoundaries(obj, onlyMaterialInterfaces)
            if ~exist("onlyMaterialInterfaces", "var")
                onlyMaterialInterfaces = false;
            end

            LWnomat = 0.5;
            if onlyMaterialInterfaces
                LWmat = 1.5;
            else
                LWmat = 0.5;
            end
            
            obj.checkRotationUpdate();

            % plot certain interfaces
            for iInt = 1:numel(obj.Interfaces)
                % Check if interface is on excitation boundary
                excitationFlag = false;
                for iEXC = 1:numel(obj.Excitations)
                    if xor(any(obj.Excitations(iEXC).Excitation.Patches == obj.Interfaces(iInt).patch1),...
                            any(obj.Excitations(iEXC).Excitation.Patches == obj.Interfaces(iInt).patch2))
                        excitationFlag = true;
                        break;
                    end
                    % Check for different phases
                    if isa(obj.Excitations(iEXC).Excitation, "IGA_EXCcurrent")
                        exc = obj.Excitations(iEXC).Excitation;
                        for iPhase = 1:numel(exc.Phases)
                            if xor(any(exc.Phases(iPhase).Patches == obj.Interfaces(iInt).patch1),...
                                    any(exc.Phases(iPhase).Patches == obj.Interfaces(iInt).patch2))
                                excitationFlag = true;
                                break;
                            end
                        end
                    end
                end
                % Select linewidth depending if material boundaries should be drawn
                if obj.Geometry(obj.Interfaces(iInt).patch1).Material ~= obj.Geometry(obj.Interfaces(iInt).patch2).Material
                    isMaterialInterface = true;
                elseif excitationFlag
                    isMaterialInterface = true;
                else
                    isMaterialInterface = false;
                end
                if isMaterialInterface
                    LW = LWmat;
                    COL = [0,0,0]; % obj.Interfaces(iInt).PlotColor;
                else
                    if onlyMaterialInterfaces
                        continue
                    end
                    LW = LWnomat;
                    COL = obj.Interfaces(iInt).PlotColor;
                end
                if ishandle(obj.FigureGeometryInterfaces(iInt))
                    set(obj.FigureGeometryInterfaces(iInt), 'XData', obj.Interfaces(iInt).PlotPoints{1}, 'YData', obj.Interfaces(iInt).PlotPoints{2}, 'Color', COL, 'LineWidth', LW);
                else
                    obj.FigureGeometryInterfaces(iInt) = plot(obj.Interfaces(iInt).PlotPoints{1}, obj.Interfaces(iInt).PlotPoints{2}, "Color", COL, "LineWidth", LW);
                end
            end
            % plot all boundaries
            for iBnd = 1:numel(obj.Boundaries)
                if ishandle(obj.FigureGeometryBoundaries(iBnd))
                    set(obj.FigureGeometryBoundaries(iBnd), 'XData', obj.Boundaries(iBnd).PlotPoints{1}, 'YData', obj.Boundaries(iBnd).PlotPoints{2});
                else
                    obj.FigureGeometryBoundaries(iBnd) = plot(obj.Boundaries(iBnd).PlotPoints{1}, obj.Boundaries(iBnd).PlotPoints{2}, "Color", obj.Boundaries(iBnd).PlotColor, "LineWidth", LWmat);
                end
            end
        end

        function exportIGES(obj, filename, scaling)
            if ~exist("filename", "var")
                filename = 'nrb.igs';
            end
            if ~exist("scaling", "var")
                scaling = 1;
            else
                disp("SCALING NOT YET INCLUDED")
            end

            % Add all NURBS boundaries to export
            for iBnd = 1:numel(obj.Boundaries)
                patch = obj.Boundaries(iBnd).patches;
                side = obj.Boundaries(iBnd).faces;
                nrb_export(iBnd) = nrbextract(obj.Geometry(patch).nurbs, side);
            end

            % Add NURBS interfaces, which have different materials or excitations
            for iInt = 1:numel(obj.Interfaces)
                % Check if interface is on excitation boundary
                excitationFlag = false;
                for iEXC = 1:numel(obj.Excitations)
                    if xor(any(obj.Excitations(iEXC).Excitation.Patches == obj.Interfaces(iInt).patch1),...
                            any(obj.Excitations(iEXC).Excitation.Patches == obj.Interfaces(iInt).patch2))
                        excitationFlag = true;
                        break;
                    end
                    % Check for different phases
                    if isa(obj.Excitations(iEXC).Excitation, "IGA_EXCcurrent")
                        exc = obj.Excitations(iEXC).Excitation;
                        for iPhase = 1:numel(exc.Phases)
                            if xor(any(exc.Phases(iPhase).Patches == obj.Interfaces(iInt).patch1),...
                                    any(exc.Phases(iPhase).Patches == obj.Interfaces(iInt).patch2))
                                excitationFlag = true;
                                break;
                            end
                        end
                    end
                end
                % Select linewidth depending if material boundaries should be drawn
                if obj.Geometry(obj.Interfaces(iInt).patch1).Material ~= obj.Geometry(obj.Interfaces(iInt).patch2).Material
                    isMaterialInterface = true;
                elseif excitationFlag
                    isMaterialInterface = true;
                else
                    isMaterialInterface = false;
                end
                if isMaterialInterface
                    patch = obj.Interfaces(iInt).patch1;
                    side = obj.Interfaces(iInt).side1;
                    nrb_export(end+1) = nrbextract(obj.Geometry(patch).nurbs, side);
                end
            end
            % Export to IGES
            scaling = 1000;
            for inrb = 1:numel(nrb_export)
                nrb_export(inrb).coefs(1:2,:) = nrb_export(inrb).coefs(1:2,:)*scaling;
            end
            nrb2iges(nrb_export, filename)
        end

        function plotBoundaryConditions(obj)
            lw = 2;
            for iBC = 1:numel(obj.BoundaryConditions)
                checkBC = obj.BoundaryConditions(iBC).BoundaryCondition;
                if isprop(checkBC, "BoundaryNumbers")
                    col = checkBC.getPlotColor();
                    for bndNr =  checkBC.BoundaryNumbers
                        plot(obj.Boundaries(bndNr).PlotPoints{1}, obj.Boundaries(bndNr).PlotPoints{2}, "Color", col, "LineWidth", lw);
                    end
                end
            end
        end

        function plotPatchNumbers(obj)
            for iPatch = 1:obj.NumberPatches
                x = obj.Meshes.msh_patch{iPatch}.map([0.5, 0.5]);
                text(x(1), x(2), num2str(iPatch), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
            end
        end

        function plotBoundaryNumbers(obj)
            for iBnd = 1:numel(obj.Boundaries)
                x = nrbeval(obj.Geometry(obj.Boundaries(iBnd).patches).boundary(obj.Boundaries(iBnd).faces).nurbs, 0.5);
                text(x(1), x(2), num2str(iBnd), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
            end
        end

        function plotInterfaceNumbers(obj)
            for iInt = 1:numel(obj.Interfaces)
                x = nrbeval(obj.Geometry(obj.Interfaces(iInt).patch1).boundary(obj.Interfaces(iInt).side1).nurbs, 0.5);
                text(x(1), x(2), num2str(iInt), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
            end
        end

        function plotBoundaryInterfaceSides(obj, bndNr)
            patch = obj.Boundaries(bndNr).patches;
            side = obj.Boundaries(bndNr).faces;
            x1 = nrbeval(obj.Geometry(patch).boundary(side).nurbs, 0);
            x2 = nrbeval(obj.Geometry(patch).boundary(side).nurbs, 1);
            text(x1(1), x1(2), "1", "HorizontalAlignment", "center", "VerticalAlignment", "middle");
            text(x2(1), x2(2), "2", "HorizontalAlignment", "center", "VerticalAlignment", "middle");
        end

        function moveInterfaceToBoundary(obj, interfaceNr)
            for iInt = 1:numel(interfaceNr)
                obj.Boundaries(end+1).nsides = 1;
                obj.Boundaries(end).patches = obj.Interfaces(interfaceNr(iInt)).patch1;
                obj.Boundaries(end).faces = obj.Interfaces(interfaceNr(iInt)).side1;

                obj.Boundaries(end+1).nsides = 1;
                obj.Boundaries(end).patches = obj.Interfaces(interfaceNr(iInt)).patch2;
                obj.Boundaries(end).faces = obj.Interfaces(interfaceNr(iInt)).side2;
            end
            obj.Interfaces(interfaceNr) = [];
            warning("BoundaryInterfaces might not be correct, you need to add them manually!")
        end

        function plotControlPoints(obj, drawText, indices, markerSize, col)
            if ~exist("drawText", "var") || isempty(drawText)
                drawText = false;
            end
            if ~exist("indices", "var")  || isempty(indices)
                indices = 1:size(obj.ControlPoints, 1);
            end
            if ~exist("markerSize", "var") || isempty(markerSize)
                markerSize = 15;
            end
            if ~exist("col", "var") || isempty(col)
                col = TUDa_getColor("9b");
            end

            scatter(obj.ControlPoints(indices, 1), obj.ControlPoints(indices, 2), "filled", 'MarkerFaceColor', col, 'SizeData', markerSize, Marker='o')
            if drawText
                text(obj.ControlPoints(indices, 1), obj.ControlPoints(indices, 2), string(indices), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
            end
        end

        function plotQuadraturePoints(obj, patches, markerSize, col)
            if ~exist("patches", "var") || isempty(patches)
                patches = 1:obj.NumberPatches;
            end
            if ~exist("markerSize", "var") || isempty(markerSize)
                markerSize = 15;
            end
            if ~exist("col", "var")
                col = TUDa_getColor("9b");
            end
            points = [];
            for patch = patches
                points = [points, reshape(obj.MeshesEval{patch}.geo_map, 2, [])];

            end
            scatter(points(1,:), points(2,:), "filled", 'MarkerFaceColor', col, 'SizeData', markerSize, Marker='o');
        end

        function plotEvaluationPoints(obj, plotNumbers)
            if ~exist("plotNumbers", "var")
                plotNumbers = false;
            end
            xPhysical = zeros(size(obj.EvaluationPoints));
            for iPoint = 1:size(obj.EvaluationPoints, 2)
                xPhysical(:,iPoint) = obj.Meshes.msh_patch{obj.EvaluationPatches(iPoint)}.map(obj.EvaluationPoints(:, iPoint));
            end
            scatter(xPhysical(1,:), xPhysical(2,:), "black", "filled");
            if plotNumbers
                text(xPhysical(1,:), xPhysical(2,:), " " +string(1:size(xPhysical,2)), "HorizontalAlignment", "left", "VerticalAlignment", "middle");
            end
        end

        % Calculate the Shape Values and Gradients at Plot Points in advance
        function calculatePlotValues(obj)
            progressbar([obj.Name ': Calculating function space at plot points: '])
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).PlotValues = shp_eval1(obj.Spaces.sp_patch{iPatch}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsParam , {'value', 'gradient'});
                progressbar(iPatch/obj.NumberPatches);
            end
            progressbar('finished');
        end

        function initializeFigures(obj)
            obj.FigureGeometry = gobjects(1);
            obj.FigureGeometryInterfaces = gobjects(numel(obj.Interfaces), 1);
            obj.FigureGeometryBoundaries = gobjects(numel(obj.Boundaries), 1);
            obj.FigureGeometryPatches = gobjects(obj.NumberPatches, 1);
        end

        function pnts = getGeoPlotPoints(~, nurbs, pts, rot)
            if nargin <= 3
                rot = 0;
            end
            xpts = pts{1};
            ypts = pts{2};

            F_geo = nrbeval(nurbs, {xpts, ypts});
            X = squeeze(F_geo(1, :, :))*cos(rot) - squeeze(F_geo(2, :, :))*sin(rot);
            Y = squeeze(F_geo(2, :, :))*cos(rot) + squeeze(F_geo(1, :, :))*sin(rot);
            pnts{1, 1} = X;
            pnts{2, 1} = Y;

            bnd = nrbeval(nrbextract(nurbs, 1), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 2} = X;
            pnts{2, 2} = Y;

            bnd = nrbeval(nrbextract(nurbs, 2), ypts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 3} = X;
            pnts{2, 3} = Y;

            bnd = nrbeval(nrbextract(nurbs, 3), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 4} = X;
            pnts{2, 4} = Y;

            bnd = nrbeval(nrbextract(nurbs, 4), xpts);
            X = bnd(1, :)*cos(rot) - bnd(2, :)*sin(rot);
            Y = bnd(2, :)*cos(rot) + bnd(1, :)*sin(rot);
            pnts{1, 5} = X;
            pnts{2, 5} = Y;
        end

        function exportXML(obj, name)
            disp("THIS IS WORK IN PROGRESS: connectivity and nurbs weights not yet included")
            fileID = fopen(name, 'w');
            fprintf(fileID, '<xml>\n');
            arr2text = @(text) regexprep(num2str(text),'\s+',' ');
            for iPatch = 1:obj.NumberPatches
                fprintf(fileID, '  <!-- Patch %i -->\n', iPatch);
                type = 'TensorBSpline2';
                id = iPatch;
                basis2d = 'TensorBSplineBasis2';
                basis1d = 'BSplineBasis';
                fprintf(fileID, '  <Geometry type="%s" id="%i">\n', type, id);
                fprintf(fileID, '    <Basis type="%s">\n', basis2d);
                fprintf(fileID, '      <Basis type="%s" index="0">\n', basis1d);
                fprintf(fileID, '        <KnotVector degree="%i">%s </KnotVector>\n', obj.Geometry(iPatch).nurbs.order(1)-1, arr2text(obj.Geometry(iPatch).nurbs.knots{1}));
                fprintf(fileID, '      </Basis>\n');
                fprintf(fileID, '      <Basis type="%s" index="1">\n', basis1d);
                fprintf(fileID, '        <KnotVector degree="%i">%s </KnotVector>\n', obj.Geometry(iPatch).nurbs.order(2)-1, arr2text(obj.Geometry(iPatch).nurbs.knots{2}));
                fprintf(fileID, '      </Basis>\n');
                fprintf(fileID, '    </Basis>\n');
                fprintf(fileID, '    <coefs geoDim="%i">\n', obj.Geometry(iPatch).rdim);
                coefs = reshape(obj.Geometry(iPatch).nurbs.coefs(1:obj.Geometry(iPatch).rdim,:,:), obj.Geometry(iPatch).rdim, [])';
                for iCtrlpnt = 1:size(coefs, 1)
                    fprintf(fileID, '      %s\n', arr2text(coefs(iCtrlpnt, :)));
                end
                fprintf(fileID, '    </coefs>\n');
                fprintf(fileID, '  </Geometry>\n');
                fprintf(fileID, '\n');
            end
            % interfaces
            fprintf(fileID, '  <!-- Multipatch -->\n');
            parDim = 2;
            id = 0;
            fprintf(fileID, '  <MultiPatch parDim="%i" id="%i">\n', parDim, id);
            fprintf(fileID, '    <patches type="id_range">%i %i</patches>\n', 1, obj.NumberPatches);
            fprintf(fileID, '    <interfaces>\n');
            for iInt = 1:numel(obj.Interfaces)
                interfaces = [obj.Interfaces(iInt).patch1, obj.Interfaces(iInt).side1, obj.Interfaces(iInt).patch2, obj.Interfaces(iInt).side2, [0,1,1,1]];
                fprintf(fileID, '      %s \n', arr2text(interfaces)); %1 2 0 1 0 1 1 1
            end

            fprintf(fileID, '    </interfaces>\n');

            fprintf(fileID, '    <boundary>\n');
            for iBnd = 1:numel(obj.Boundaries)
                boundary = [obj.Boundaries(iBnd).patches, obj.Boundaries(iBnd).faces];
                fprintf(fileID, '      %s \n', arr2text(boundary));
            end
            fprintf(fileID, '    </boundary>\n');
            fprintf(fileID, '  </MultiPatch>\n');

            fprintf(fileID, '</xml>');
            fclose(fileID);
        end
    
        function checkRotationUpdate(obj)
            if obj.RotationAngle ~= obj.PlotAngle
                obj.updatePlotPoints(obj.RotationAngle);
                obj.PlotAngle = obj.RotationAngle;
            end
        end
    end
end
