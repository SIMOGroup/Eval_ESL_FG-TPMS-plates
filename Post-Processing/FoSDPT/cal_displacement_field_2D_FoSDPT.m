function [FE_disp] = cal_displacement_field_2D_FoSDPT(IGA,Plate,z_coords)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate displacement field on FE grid of rectangular plate with 57-Dof 2D SS's FoSDPT plate models %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
p = IGA.NURBS.p; uKnot = IGA.NURBS.uKnot; mcp = IGA.NURBS.mcp;
q = IGA.NURBS.q; vKnot = IGA.NURBS.vKnot; ncp = IGA.NURBS.ncp;
% CP = NURBS.CP;
Ien = IGA.NURBS.Ien; Inn = IGA.NURBS.Inn;
% Iee = IGA.NURBS.Iee; Ine = IGA.NURBS.Ine;
FE_grid = IGA.NURBS.FE_grid;
ndof = IGA.params.ndof;
U = IGA.result.U;

%% Used parameters from Plate
h = Plate.geo.h;

%% ===== Initial displacement matrices =====
npoint_x = size(FE_grid,1); npoint_y = size(FE_grid,2); npoint_z = length(z_coords);
FE_disp = zeros(npoint_x,npoint_y,npoint_z,3);

%% ===== Displacement field =====
for j = 1:npoint_y
    v = vKnot(q+1) + (vKnot(ncp+1) - vKnot(q+1))/(npoint_y - 1) * (j-1);
    nj = find_knot_span(v,vKnot,ncp);
    for i = 1:npoint_x
        u = uKnot(p+1) + (uKnot(mcp+1) - uKnot(p+1))/(npoint_x - 1) * (i-1);
        ni = find_knot_span(u,uKnot,mcp);
        
        % --- Element parameters ---
        idx = find(Inn(:,1) == ni & Inn(:,2) == nj);
        c_ele = (Ien(:,1) == idx);
        sctr = Ien(c_ele,:);  % Control points indexes
        nn = length(sctr);  % Number of control points in the element
        for idof = 1:ndof  % Dofs of control points
            sctrF(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
        end
        
        % --- Kinematic matrices of NURBS shape function ---
        [N, ~, ~, dNdx, ~, ~] = cal_Kine_Shape_2D_2nd(IGA.NURBS,ni,nj,u,v);
        
        N0 = zeros(3,ndof*nn);
        N0(1,1:ndof:ndof*nn) = N';
        N0(2,2:ndof:ndof*nn) = N';
        N0(3,3:ndof:ndof*nn) = N';
        
        N1 = zeros(3,ndof*nn);
        N1(1,3:ndof:ndof*nn) = -dNdx(:,1)';
        N1(2,3:ndof:ndof*nn) = -dNdx(:,2)';
        
        N2 = zeros(3,ndof*nn);
        N2(1,4:ndof:ndof*nn) = N';
        N2(2,5:ndof:ndof*nn) = N';
        
        N3 = zeros(3,ndof*nn);
        N3(1,6:ndof:ndof*nn) = N';
        N3(2,7:ndof:ndof*nn) = N';
        
        N = {N0, N1, N2, N3}; 
        
        for i_z = 1:npoint_z
            z = z_coords(i_z);
            
            h_u = [1, 4/3*(z^3)/h^2, z - 4/3*(z^3)/h^2, z^2 + 2*(z^4)/h^2];

            % --- Calculation ---
            point_disp = zeros(3,1);
            for k = 1:4
                point_disp = point_disp + h_u(k)*N{k}*U(sctrF,1);
            end

            FE_disp(i,j,i_z,:) = point_disp;
        end
    end
end
end
