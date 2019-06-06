function [CostF,GradF,HessF,Pas1F,Pas2F,GrP1F,GrP2F] = PI_OptTools(basis_mats)
% PI_OPTTOOLS:      Construct the cost, gradient, Hessian, passivity
%                   constraint, and passivity constraint gradient functions 
%                   which are associated with the minimization function. 
%                   Add regularized prior. This function assumes that the 
%                   basis matrices *not* included in the "basis_mats" 
%                   structure are also not going to be included in the 
%                   prior matrix.
%
% Inputs:         
% 1) basis_mats   List of basis matrices which will be used in constructing
%                 the admittance matrix. There are 8 options, where 1
%                 corresponds to inclusion and 0 corresponds to exclusion
%                     1) basis_mats.T1
%                     2) basis_mats.T2
%                     3) basis_mats.T3
%                     4) basis_mats.T4
%                     5) basis_mats.jT1
%                     6) basis_mats.jT2
%                     7) basis_mats.jT3
%                     8) basis_mats.jT4
%
% Outputs:
% 1) CostF        Cost function
% 2) GradF        Gradient function
% 3) HessF        Hessian function
% 4) Pas1F        Passivity constraint (eigenvalue 1)
% 5) Pas2F        Passivity constraint (eigenvalue 2)
% 6) GrP1F        Passivity constraint gradient (eigenvalue 1)
% 7) GrP2F        Passivity constraint gradient (eigenvalue 2)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = sqrt(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note - The objective function is scaled by 1e5 for numerical reasons %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define basis matrices
T1 = [1  0; 0  1];
T2 = [1  0; 0 -1];
T3 = [0  1; 1  0];
T4 = [0 -1; 1  0];

% Define parameter vector
param_vec = [];
prior_vec = [];

% Define symbolic coefficient variables if applicable: T1
if basis_mats.T1 == 1
    aT1       = sym('aT1','real');
    aT1p      = sym('aT1p','real');
    param_vec = [param_vec; aT1];
    prior_vec = [prior_vec; aT1p];
else
    aT1  = 0;
    aT1p = 0;
end

% T2
if basis_mats.T2 == 1
    aT2       = sym('aT2','real');
    aT2p      = sym('aT2p','real');
    param_vec = [param_vec; aT2];
    prior_vec = [prior_vec; aT2p];
else
    aT2  = 0;
    aT2p = 0;
end

% T3
if basis_mats.T3 == 1
    aT3       = sym('aT3','real');
    aT3p      = sym('aT3p','real');
    param_vec = [param_vec; aT3];
    prior_vec = [prior_vec; aT3p];
else
    aT3  = 0;
    aT3p = 0;
end

% T4
if basis_mats.T4 == 1
    aT4       = sym('aT4','real');
    aT4p      = sym('aT4p','real');
    param_vec = [param_vec; aT4];
    prior_vec = [prior_vec; aT4p];
else
    aT4  = 0;
    aT4p = 0;
end

% jT1
if basis_mats.jT1 == 1
    bT1       = sym('bT1','real');
    bT1p      = sym('bT1p','real');
    param_vec = [param_vec; bT1];
    prior_vec = [prior_vec; bT1p];
else
    bT1  = 0;
    bT1p = 0;
end

% jT2
if basis_mats.jT2 == 1
    bT2       = sym('bT2','real');
    bT2p      = sym('bT2p','real');
    param_vec = [param_vec; bT2];
    prior_vec = [prior_vec; bT2p];
else
    bT2  = 0;
    bT2p = 0;
end

% jT3
if basis_mats.jT3 == 1
    bT3       = sym('bT3','real');
    bT3p      = sym('bT3p','real');
    param_vec = [param_vec; bT3];
    prior_vec = [prior_vec; bT3p];
else
    bT3  = 0;
    bT3p = 0;
end

% jT4
if basis_mats.jT4 == 1
    bT4       = sym('bT4','real');
    bT4p      = sym('bT4p','real');
    param_vec = [param_vec; bT4];
    prior_vec = [prior_vec; bT4p];
else
    bT4  = 0;
    bT4p = 0;
end

% Build admittance matrix Y and parse
Y   = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4 + j*bT1*T1 + j*bT2*T2 + j*bT3*T3 + j*bT4*T4;
Y11 = Y(1,1);
Y12 = Y(1,2);
Y21 = Y(2,1);
Y22 = Y(2,2);

% Define complex symbolic inputs and outputs
syms y1 y2 u1 u2

% Regularizer
lambda = sym('lambda','real');

% Prior Function
fp = lambda*((aT1-aT1p)^2 + (aT2-aT2p)^2 + (aT3-aT3p)^2 + (aT4-aT4p)^2 + ...
             (bT1-bT1p)^2 + (bT2-bT2p)^2 + (bT3-bT3p)^2 + (bT4-bT4p)^2);

% Cost Function
f = (y1 - Y11*u1 - Y12*u2)*conj(y1 - Y11*u1 - Y12*u2) + ...
    (y2 - Y21*u1 - Y22*u2)*conj(y2 - Y21*u1 - Y22*u2);

% Objective Function
obj_func = 1e4*(f+fp);

% Take symbolic gradient and Hessian: take real part to ensure
% numerically small imaginary part does not interfere
CostS = real(obj_func);
GradS = real(jacobian(obj_func,param_vec));
HessS = real(hessian(obj_func,param_vec));

% Passivity Constraints
Pas1S = (bT4 + sqrt(aT1^2 + bT2^2 + bT3^2));
Pas2S = (bT4 - sqrt(aT1^2 + bT2^2 + bT3^2));

% Passivity Gradients
GrP1S = jacobian(Pas1S,param_vec);
GrP2S = jacobian(Pas2S,param_vec);

% Concatenate Variables
input  = [u1; u2];
output = [y1; y2];

% Transform into matlab functions
CostF = matlabFunction(CostS,'vars',{input,output,lambda,prior_vec,param_vec});
GradF = matlabFunction(GradS,'vars',{input,output,lambda,prior_vec,param_vec});
HessF = matlabFunction(HessS,'vars',{input,output,lambda,prior_vec,param_vec});
Pas1F = matlabFunction(Pas1S,'vars',{param_vec});
Pas2F = matlabFunction(Pas2S,'vars',{param_vec});
GrP1F = matlabFunction(GrP1S,'vars',{param_vec});
GrP2F = matlabFunction(GrP2S,'vars',{param_vec});

end

