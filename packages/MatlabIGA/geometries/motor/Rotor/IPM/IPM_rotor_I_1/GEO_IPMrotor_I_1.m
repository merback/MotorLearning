% https://www.mathcha.io/editor/P1GLZsWPU5OtxQ7E6rh1q5B6PTDw7LEpTExPn4G

function [srf, patches, params] = GEO_IPMrotor_I_1(parameters)
    % Standard parameters
    params.draw_geometry = true;
    params.POLES       = 6;
    params.Symmetry    = 6;                 % Reduction due to symmetry
    params.MW          = 19e-3;             % Magnet Width
    params.MH          = 7e-3;              % Magnet Height
    params.MAG         = 7e-3;              % Magnet to Air Gap distance
    params.RRO         = 44e-3;             % Rotor Radius Outer
    params.RRI         = 16e-3;             % Rotor Radius Inner
    params.RAG         = 44.5e-3;           % Air Gap Radius
    params.Theta1      = 0.37;              % slit positioning
    params.Theta2      = 0.37558;           % slit positioning
    params.LW1         = 2e-3;              % slit positioning
    params.LW2         = 3e-3;              % slit positioning
    params.REFINEMENT  = 7;
    params.Type        = "nonlinear";
    % Overwrite standard parameters with input parameters
    if nargin == 1
        param_names = fieldnames (parameters);
        for iParam  = 1:numel (param_names)
            if isfield(params, param_names(iParam))
                eval(['params.' param_names{iParam} '= parameters.(param_names{iParam});']);
            end
        end
    end

    if params.Type == "linear"
        patches.Iron = [1,3,4,6,9,10,11];
        patches.Magnets = [8];
        patches.Air = [7,12:17,2,5];
    elseif params.Type == "nonlinear"
        patches.Iron = [1,2,3,4,5,6,9,10,11];
        patches.Magnets = [8];
        patches.Air = [7,12:17];
    else
        error("Wrong Type");
    end
    
    % Derived parameters
    R_outRt = params.RRO;
    R_inRt = params.RRI;
    R_gamma = params.RAG;
    beta = 2*pi/params.POLES;
    alpha = beta/2;
    center = [0, 0, 0];
    dMag = R_outRt - params.MH - params.MAG;
    r6 = R_outRt - params.LW2;
    r7 = R_outRt - params.LW1;
    THETA1 = 0.35;
    THETA2 = 0.4;
    
    % setting points
    point{1} = R_inRt*[cos(alpha), -sin(alpha)];
    point{2} = R_inRt*[cos(alpha), sin(alpha)];
    
    point{3} = [dMag, -params.MW/2];
    point{4} = [dMag, params.MW/2];
    
    point{5} = r6*[cos(alpha), -sin(alpha)];
    point{6} = r6*[cos(params.Theta2), -sin(params.Theta2)];
    point{7} = r7*[cos(params.Theta1), -sin(params.Theta1)];

    point{8} = [dMag + params.MH, -params.MW/2];
    point{9} = [dMag + params.MH, params.MW/2];

    point{10} = r7*[cos(params.Theta1), sin(params.Theta1)];    
    point{11} = r6*[cos(params.Theta2), -sin(params.Theta2)];
    point{12} = r6*[cos(alpha), -sin(alpha)];

    point{13} = R_outRt*[cos(alpha), -sin(alpha)];
    point{14} = R_outRt*[cos(params.Theta2), -sin(params.Theta2)];
    point{15} = R_outRt*[cos(params.Theta1), -sin(params.Theta1)];
    point{16} = R_outRt*[cos(params.Theta1), sin(params.Theta1)];
    point{17} = R_outRt*[cos(params.Theta2), sin(params.Theta2)];
    point{18} = R_outRt*[cos(alpha), sin(alpha)];
    
    point{19} = R_gamma*[cos(alpha), -sin(alpha)];
    point{20} = R_gamma*[cos(THETA2), -sin(THETA2)];
    point{21} = R_gamma*[cos(THETA1), -sin(THETA1)];
    point{22} = R_gamma*[cos(THETA1), sin(THETA1)];
    point{23} = R_gamma*[cos(THETA2), sin(THETA2)];
    point{24} = R_gamma*[cos(alpha), sin(alpha)];
    
    % Patch 1 - Bottom Right Iron Part
    LeftSide = nrbline(point{5}, point{6});
    RightSide  = nrbcirc (R_outRt, center, -alpha, -params.Theta2);
    srf(1) = nrbruled(LeftSide, RightSide);

    % Patch 2 - Bottom Middle Air Slit
    LeftSide = nrbline(point{6}, point{7});
    TopSide = nrbline(point{7}, point{15});
    RightSide = nrbcirc(R_outRt, center, -params.Theta2, -params.Theta1);
    BotSide = nrbline(point{6}, point{14});
    srf(2) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 3 - Right Iron Part
    LeftSide = nrbline(point{7}, point{10});
    RightSide  = nrbcirc (R_outRt, center, -params.Theta1, params.Theta1);
    refinedKnot = kntrefine(RightSide.knots, params.REFINEMENT, RightSide.order-1, RightSide.order-2);
    RightSide = nrbkntins(RightSide, setdiff(refinedKnot, RightSide.knots));

    srf(3) = nrbruled(LeftSide, RightSide);

    % Patch 4 - Middle Iron Part
    LeftSide = nrbline(point{8}, point{9});
    TopSide = nrbline(point{9}, point{10});
    RightSide = nrbextract(srf(3), 3); % = nrbline(point{7}, point{10});
    BotSide = nrbline(point{8}, point{7});
    srf(4) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    
    mirrorx = diag([1,-1,1,1]);
    % Patch 5 - Top Middle Air Slit
    srf(5) = nrbtform(srf(2), mirrorx);

    % Patch 6 - Top Right Iron Part
    srf(6) = nrbtform(srf(1), mirrorx);

    % Patch 7 - Bottom Left Air Slit
    LeftSide = nrbline(point{6}, point{3});
    TopSide = nrbline(point{3}, point{8});
    RightSide = nrbreverse(nrbextract(srf(4), 3)); %nrbline(point{7}, point{8});
    BotSide = nrbextract(srf(2), 1);  %nrbline(point{6}, point{7});
    srf(7) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 8 - Permanent Magnet
    LeftSide = nrbline(point{3}, point{4});
    TopSide = nrbline(point{4}, point{9});
    RightSide = nrbextract(srf(4), 1); %nrbline(point{8}, point{9});
    BotSide = nrbextract(srf(7), 4); %nrbline(point{3}, point{8});
    srf(8) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 9 - Bottom Left Iron Part
    LeftSide = nrbline(point{1}, point{3});
    TopSide = nrbline(point{3}, point{6});
    RightSide = nrbextract(srf(1), 3); %nrbline(point{5}, point{6});
    BotSide = nrbline(point{1}, point{5});
    srf(9) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    
    % Patch 10 - Left Iron Part
    TopSide = nrbline(point{2}, point{4});
    LeftSide  = nrbcirc(R_inRt, center, -alpha, alpha);
    RightSide = nrbextract(srf(8), 1); %nrbline(point{3}, point{4});
    BotSide = nrbextract(srf(9), 1);
    srf(10) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 11 - Top Left Iron Part
    srf(11) = nrbtform (srf(9), mirrorx);
    
    % Patch 12 - Top Left Air Slit
    srf(12) = nrbtform(srf(7), mirrorx);

    % Patch 13 - Bottom Right Air Slit
    LeftSide = nrbextract(srf(1), 4); %= nrbcirc (R_outRt, center, -alpha, -params.Theta2);
    RightSide  = nrbcirc(R_gamma, center, -alpha, -THETA2);
    srf(13) = nrbruled(LeftSide, RightSide);

    % Patch 14 - Bottom/Middle Rigth Air Slit
    LeftSide = nrbextract(srf(2), 2); %  nrbcirc(R_outRt, center, -params.Theta2, -params.Theta1);
    RightSide  = nrbcirc(R_gamma, center, -THETA2, -THETA1);
    srf(14) = nrbruled(LeftSide, RightSide);

    % Patch 15 - Middle Rigth Air Slit
    LeftSide = nrbextract(srf(3), 4); % nrbcirc(R_outRt, center, -params.Theta1, params.Theta1);
    RightSide  = nrbcirc(R_gamma, center, -THETA1, THETA1);
    srf(15) = nrbruled(LeftSide, RightSide);

    % Patch 15 - Middle/Top Rigth Air Slit
    srf(16) = nrbtform(srf(14), mirrorx);

    % Patch 17 - Top Right Air Slit
    srf(17) = nrbtform(srf(13), mirrorx);

       
    % for iptc = 1:numel (srf)
    %   srf(iptc) = nrbtform (srf(iptc), vecrotz(pi/2));
    % end

    % Copy and rotate
    nsrf = numel(srf);
    IronInit = patches.Iron;
    AirInit = patches.Air;
    MagnetsInit = patches.Magnets;
    for I = 1:params.POLES/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/params.POLES));
        end
        patches.Iron = [patches.Iron, IronInit+nsrf*I];
        patches.Air = [patches.Air, AirInit+nsrf*I];
        patches.Magnets = [patches.Magnets, MagnetsInit+nsrf*I];
    end
    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for I = 1:nsrf
        srf(I) = nrbtform(srf(I), vecrotz(pi/params.POLES + pi/2 - pi/params.Symmetry));
    end

    if (params.draw_geometry)
        figure()
        clf
        hold on
        axis equal
        for i = 1:numel(srf)
            if ismember(i, patches.Iron)
                nrbplotcol(srf(i), [10,10], 'color', [0.3,0.3,0.3]);
            elseif ismember(i, patches.Magnets)
                nrbplotcol(srf(i), [10,10], 'color', 'green');
            elseif ismember(i, patches.Air)
                nrbplotcol(srf(i), [10,10], 'color', [0.1,0.1,1]);
            elseif ismember(i, patches.Windings)
                nrbplotcol(srf(i), [10,10], 'color', [0.9,0.1,1]);
            else
                disp("Rotor patch is not material defined")
                disp(i)
                nrbplotcol(srf(i), [10,10], 'color', [0.1,0.1,0.1]);
            end
        end

        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+1e-5, point{i}(2), string(i))
        % end
        view(2)
        axis equal
        drawnow();
    end
end
