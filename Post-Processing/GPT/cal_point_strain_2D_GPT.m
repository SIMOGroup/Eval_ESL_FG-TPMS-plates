function [point_strain] = cal_point_strain_2D_GPT(IGA,Plate,phys_coord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate strain field of specific point of rectangular plate with 5-Dof 2D GPT plate models %%%
%%% Strain = {\epsilon_xx, \epsilon_yy, \gamma_xy, \gamma_yz, \gamma_zx, \epsilon_zz}^T %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Ien = IGA.NURBS.Ien; Inn = IGA.NURBS.Inn;
Iee = IGA.NURBS.Iee; Ine = IGA.NURBS.Ine;
gcoord = IGA.NURBS.gcoord; nel = IGA.NURBS.nel;
ndof = IGA.params.ndof;
U = IGA.result.U;

%% Used parameters from Plate
h = Plate.geo.h;
shear_func = Plate.theory.shear_func; n_iterated = Plate.theory.iterated_order;

%% ===== Initial displacement matrices =====
point_strain = zeros(6,1);

%% ===== Point displacement =====
% --- Find element ---
for i_ele = 1:size(Iee,1)    
    % Search within element coordinate box
    edges = Iee(i_ele, :);
    box_min = reshape(min(min(Ine(edges,:,:),[],1),[],2),1,[]);
    box_max = reshape(max(max(Ine(edges,:,:),[],1),[],2),1,[]);
    
    for i_dim = 1:2
        if phys_coord(i_dim) >= box_min(i_dim) && phys_coord(i_dim) <= box_max(i_dim)
            check = true;
        else
            check = false;
            break
        end
    end
    if check
        c_ele = i_ele;
        break
    end
end

% --- Element parameters ---
sctr = Ien(c_ele,:);  % Control points indexes
ni = Inn(Ien(c_ele,1),1);  % Index of the element in parametric domain
nj = Inn(Ien(c_ele,1),2);
nn = length(sctr);  % Number of control points in the element
for idof = 1:ndof  % Dofs of control points
    sctrF(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
end
nodes = gcoord(sctr,:);

% --- Point parameters ---
nat_coord = cal_nat_coord_surf(IGA.NURBS,ni,nj,nodes,phys_coord);  % Natural coordinate
nat_coord(1) = (uKnot(ni+1)-uKnot(ni))/2*nat_coord(1) + (uKnot(ni+1)+uKnot(ni))/2;  % Map the point to parametric domain
nat_coord(2) = (vKnot(nj+1)-vKnot(nj))/2*nat_coord(2) + (vKnot(nj+1)+vKnot(nj))/2;

% --- Kinematic matrices of NURBS shape function ---
[N, ~, ~, dNdx, d2Ndx2, ~] = cal_Kine_Shape_2D_2nd(IGA.NURBS,ni,nj,nat_coord(1),nat_coord(2));

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

B2 = zeros(3,ndof*nn);  % Shear bending
B2(1,4:ndof:ndof*nn) = dNdx(:,1)';
B2(2,5:ndof:ndof*nn) = dNdx(:,2)';
B2(3,4:ndof:ndof*nn) = dNdx(:,2)';
B2(3,5:ndof:ndof*nn) = dNdx(:,1)';

Bb = {B0, B1, B2};

% Shear stiffness matrix
Bs = zeros(2,ndof*nn);  % Pure shear
Bs(1,5:ndof:ndof*nn) = N';
Bs(2,4:ndof:ndof*nn) = N';

% --- Shear deformation function ---
[fz, dfz, ~] = compute_iterated_shear_deformation_function(phys_coord(3),h,shear_func,n_iterated);
h_e = [1, phys_coord(3), fz];

% --- Calculation ---
strain_b = zeros(3,1);
for k = 1:3
    strain_b = strain_b + h_e(k)*Bb{k}*U(sctrF,1);
end
strain_s = dfz*Bs*U(sctrF,1);

point_strain = [strain_b; strain_s; 0];
point_strain = point_strain';

end
