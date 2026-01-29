function [point_disp_dist] = cal_point_displacement_distribution_2D_RPT(IGA,Plate,phys_coord,z_coords)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate displacement distribution field of specific point of rectangular plate with 4-Dof 2D RPT plate models %%%
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
npoint_z = length(z_coords);
point_disp_dist = zeros(npoint_z, 3);

%% ===== Point displacement distribution =====
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
[N, ~, ~, dNdx, ~, ~] = cal_Kine_Shape_2D_2nd(IGA.NURBS,ni,nj,nat_coord(1),nat_coord(2));

N0 = zeros(3,ndof*nn);
N0(1,1:ndof:ndof*nn) = N';
N0(2,2:ndof:ndof*nn) = N';
N0(3,3:ndof:ndof*nn) = N';
N0(3,4:ndof:ndof*nn) = N';

N1 = zeros(3,ndof*nn);
N1(1,3:ndof:ndof*nn) = -dNdx(:,1)';
N1(2,3:ndof:ndof*nn) = -dNdx(:,2)';

N2 = zeros(3,ndof*nn);
N2(1,4:ndof:ndof*nn) = dNdx(:,1)';
N2(2,4:ndof:ndof*nn) = dNdx(:,2)';

N = {N0, N1, N2}; 

for i_z = 1:npoint_z
    z = z_coords(i_z);

    % --- Shear deformation function ---
    [fz, ~, ~] = compute_iterated_shear_deformation_function(z,h,shear_func,n_iterated);
    h_u = [1, z, fz-z];
    
    % --- Calculation ---
    point_disp = zeros(3,1);
    for k = 1:3
        point_disp = point_disp + h_u(k)*N{k}*U(sctrF,1);
    end

    point_disp_dist(i_z,:) = point_disp;
end
end
