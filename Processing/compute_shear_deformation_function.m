function [fz, dfz, ddfz] = compute_shear_deformation_function(z,h,shear_func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute shear deformation function at point z %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% ["Poly"]: Polynomial, ["Trig"]: Trigonometric, ["Hype"]: Hyperbolic, ["Expo"]: Exponential, ["Comb"]: Combined or Hybrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
parts = strsplit(shear_func, '_');
func_type = parts{1}; func_idx = eval(parts{2});

%% === Shear deformation function ===
try
    switch func_type
        case "Poly"  % Polynomial
            switch func_idx
                case 1 % Reddy et al. (1984)
                    fz   =   h * ((z/h) - 4/3*(z/h)^3);
                    dfz  =   1 * (1 - 4*(z/h)^2);
                    ddfz = 1/h * (- 4*2*(z/h));
                case 2 % Nguyen-Xuan et al. (2013)
                    fz   =   h * (7/8*(z/h) - 2*(z/h)^3 + 2*(z/h)^5);
                    dfz  =   1 * (7/8 - 6*(z/h)^2 + 10*(z/h)^4);
                    ddfz = 1/h * (- 6*2*(z/h) + 10*4*(z/h)^3);
            end
        case "Trig"  % Trigonometric
            switch func_idx
                case 1 % Touratier et al. (1991)
                    fz   =   h * (1/pi*sin(pi*(z/h)));
                    dfz  =   1 * (cos(pi*(z/h)));
                    ddfz = 1/h * (-pi*sin(pi*(z/h)));
                case 2 % Thai et al. (2014)
                    fz   =   h * (atan(2*z/h) - z/h);
                    dfz  =   1 * ((1 - 4*(z/h)^2) / (1 + 4*(z/h)^2));
                    ddfz = 1/h * (-16*(z/h) / (1 + 4*(z/h)^2)^2);
                case 3 % Modified from Thai et al. (2014)
                    fz   =     1 * (atan(sin(pi*(z/h))));
                    dfz  =   1/h * (pi*cos(pi*(z/h)) / (1+(sin(pi*(z/h)))^2));
                    ddfz = 1/h^2 * (- (pi^2*sin(pi*(z/h)))/(sin(pi*(z/h))^2 + 1) - (2*pi^2*cos(pi*(z/h))^2*sin(pi*(z/h)))/((sin(pi*(z/h))^2 + 1)^2));
    
                    fz = fz*h/pi; dfz = dfz*h/pi; ddfz = ddfz*h/pi;  % Modification for function's unit
            end
        case "Hype"  % Hyperbolic
            switch func_idx
                case 1 % Meiche et al. (2011)
                    fz   =   h * ((sinh(pi*(z/h))/pi - z/h) / (cosh(pi/2) - 1) - z/h);
                    dfz  =   1 * ((cosh(pi*(z/h)) - 1) / (cosh(pi/2) - 1) - 1);
                    ddfz = 1/h * pi*sinh(pi*(z/h)) / (cosh(pi/2) - 1);
                case 2 % Modified from Thai et al. (2014)
                    fz   =     1 * (asinh(sin(pi*z/h)));
                    dfz  =   1/h * (pi*cos(pi*z/h) / sqrt(1 + sin(pi*z/h)^2));
                    ddfz = 1/h^2 * (- 2*pi^2*sin(pi*z/h) / (sin(pi*z/h)^2 + 1)^(3/2));
    
                    fz = fz*h/pi; dfz = dfz*h/pi; ddfz = ddfz*h/pi;  % Modification for function's unit
            end
        case "Expo"  % Exponential
            switch func_idx
                case 1 % Karama et al. (2003)
                    fz   =   h * z/h * exp(-2*(z/h)^2);
                    dfz  =   1 * (1 - 4*(z/h)^2) * exp(-2*(z/h)^2);
                    ddfz = 1/h * (16*(z/h)^3 - 12*(z/h)) * exp(-2*(z/h)^2);
                case 2 % Belkhodja et al. (2020)
                    fz   = pi*h / (pi^4 + h^4) * (exp(h*z/pi)*(pi^2*sin(pi*z/h) + h^2*cos(pi*z/h)) - h^2);
                    dfz  = exp(h*z/pi) * cos(pi*z/h);
                    ddfz = exp(h*z/pi) * (h/pi*cos(pi*z/h) - pi/h*sin(pi*z/h));
            end
        case "Comb"  % Combined or Hybrid
            switch func_idx
                case 1 % Modified from Mantari et al. (2012)
                    fz   =   1 * (sin(pi*(z/h)) * exp(1/2*cos(pi*z/h)) + pi*z/(2*h));
                    dfz  =   1 * (exp(cos(pi*z/h)/2) * pi/(2*h) * (cos(pi*(z/h))^2 + 2*cos(pi*(z/h)) - 1) + pi/(2*h));
                    ddfz =   1 * (- pi^2/(4*h^2) * sin(pi*z/h) * exp(cos(pi*z/h)/2) * (cos(pi*z/h)^2 + cos(pi*z/h)*6 + 3));
                    
                    fz = fz*h/(2*pi); dfz = dfz*h/(2*pi); ddfz = ddfz*h/(2*pi);  % Modification for function's unit
                case 2 % Tran et al. (2025)
                    fz   =   h * (z/h + 1/(2*pi)*sin(2*pi*(z/h)) - atan(sin(pi*(z/h))));
                    dfz  =   1 * (1 + cos(2*pi*(z/h)) - pi*cos(pi*(z/h))/(1+(sin(pi*(z/h)))^2));
                    ddfz = 1/h * (- 2*pi*sin(2*pi*z/h) - pi^2 * sin(pi*z/h) * (sin(pi*z/h)^2 - 3) / (sin(pi*z/h)^2 + 1)^2);
            end
    end
catch 
    disp('> Cannot find function')
    rethrow(ME)
end
