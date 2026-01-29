function s = find_knot_span(u,uKnot,mcp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find knot span of a natural knot point in 1D Knot vector %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source: Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: Berlin 1995
if u < uKnot(1) || u > uKnot(end)
    error('Parametric point u is not inside the Knot vector!')
end

if abs(u - uKnot(mcp+1)) < 1e-8
    s = mcp;
else
    for s = 1:length(uKnot)-1
        if u < uKnot(s+1)
            return
        end
    end  
end
end