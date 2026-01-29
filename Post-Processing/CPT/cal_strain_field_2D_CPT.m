function [FE_strain] = cal_strain_field_2D_CPT(IGA,Plate,z_coords)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate strain field on FE grid of rectangular plate with 3-Dof 2D CPT plate models %%%
%%% Strain = {\epsilon_xx, \epsilon_yy, \gamma_xy, \gamma_yz, \gamma_zx, \epsilon_zz}^T %%%
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
FE_strain = zeros(npoint_x,npoint_y,npoint_z,6);

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
        [N, ~, ~, dNdx, d2Ndx2, ~] = cal_Kine_Shape_2D_2nd(IGA.NURBS,ni,nj,u,v);
        
        % Bending stiffness matrix
        B0 = zeros(3,ndof*nn);  % Pure membrane
        B0(1,1:ndof:ndof*nn) = dNdx(:,1)';
        B0(2,2:ndof:ndof*nn) = dNdx(:,2)';
        B0(3,1:ndof:ndof*nn) = dNdx(:,2)';
        B0(3,2:ndof:ndof*nn) = dNdx(:,1)';
        
        B1 = zeros(3,ndof*nn);  % Pure bending
        B1(1,3:ndof:ndof*nn) = -d2Ndx2(:,1)';
        B1(2,3:ndof:ndof*nn) = -d2Ndx2(:,2)';
        B1(3,3:ndof:ndof*nn) = -2*d2Ndx2(:,3)';
      
        Bb = {B0, B1};
        
        for i_z = 1:npoint_z
            z = z_coords(i_z);

            h_e = [1, z];
            
            % --- Calculation ---
            strain_b = zeros(3,1);
            for k = 1:2
                strain_b = strain_b + h_e(k)*Bb{k}*U(sctrF,1);
            end
            strain_s = zeros(2,1);
                
            FE_strain(i,j,i_z,:) = [strain_b; strain_s; 0];
        end
    end
end
end
