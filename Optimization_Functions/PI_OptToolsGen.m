function [CostF,GradF,HessF,Pas1F,Pas2F,GrP1F,GrP2F] = PI_OptToolsGen(basis_mats)
% PI_OPTTOOLSGEN:   Construct the cost, gradient, Hessian, passivity
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

% Define symbolic coefficient variables if applicable: T1
if basis_mats.T1 == 1
    aT1       = sym('aT1','real');
    param_vec = [param_vec; aT1];
else
    aT1  = 0;
end

% T2
if basis_mats.T2 == 1
    aT2       = sym('aT2','real');
    param_vec = [param_vec; aT2];
else
    aT2  = 0;
end

% T3
if basis_mats.T3 == 1
    aT3       = sym('aT3','real');
    param_vec = [param_vec; aT3];
else
    aT3  = 0;
end

% T4
if basis_mats.T4 == 1
    aT4       = sym('aT4','real');
    param_vec = [param_vec; aT4];
else
    aT4  = 0;
end

% jT1
if basis_mats.jT1 == 1
    bT1       = sym('bT1','real');
    param_vec = [param_vec; bT1];
else
    bT1  = 0;
end

% jT2
if basis_mats.jT2 == 1
    bT2       = sym('bT2','real');
    param_vec = [param_vec; bT2];
else
    bT2  = 0;
end

% jT3
if basis_mats.jT3 == 1
    bT3       = sym('bT3','real');
    param_vec = [param_vec; bT3];
else
    bT3  = 0;
end

% jT4
if basis_mats.jT4 == 1
    bT4       = sym('bT4','real');
    param_vec = [param_vec; bT4];
else
    bT4  = 0;
end

% Build admittance matrix Y and parse
Y   = aT1*T1 + aT2*T2 + aT3*T3 + aT4*T4 + j*bT1*T1 + j*bT2*T2 + j*bT3*T3 + j*bT4*T4;
Y11 = Y(1,1);
Y12 = Y(1,2);
Y21 = Y(2,1);
Y22 = Y(2,2);

% Define complex symbolic inputs and outputs: source frequnecy
syms ys1 ys2 us1 us2

% Define complex symbolic inputs and outputs: non-source frequnecy
syms yns1 yns2 uns1 uns2

% Regularizer
kappa = sym('kappa','real');

% Cost Function - Source Frequnecy
f = (ys1 - Y11*us1 - Y12*us2)'*(ys1 - Y11*us1 - Y12*us2) + ...
    (ys2 - Y21*us1 - Y22*us2)'*(ys2 - Y21*us1 - Y22*us2);

% Cost Function - Non-Source Frequnecy (Regularizer Function)
fr = ((yns1 - Y11*uns1 - Y12*uns2)'*(yns1 - Y11*uns1 - Y12*uns2) + ...
      (yns2 - Y21*uns1 - Y22*uns2)'*(yns2 - Y21*uns1 - Y22*uns2));

% Objective Function: take real part to ensure numerically small 
% imaginary part does not interfere
obj_func = real(f+kappa*fr);

% Take symbolic gradient and Hessian: take real part to ensure
% numerically small imaginary part does not interfere
CostS = obj_func;
GradS = jacobian(obj_func,param_vec);
HessS = hessian(obj_func,param_vec);

% Passivity Constraints
Pas1S = (bT4 + sqrt(aT1^2 + bT2^2 + bT3^2));
Pas2S = (bT4 - sqrt(aT1^2 + bT2^2 + bT3^2));

% Passivity Gradients
GrP1S = jacobian(Pas1S,param_vec);
GrP2S = jacobian(Pas2S,param_vec);

% Concatenate Variables
input_s   = [us1; us2];
output_s  = [ys1; ys2];
input_ns  = [uns1; uns2];
output_ns = [yns1; yns2];

% Transform into matlab functions
CostF = matlabFunction(CostS,'vars',{input_s,output_s,input_ns,output_ns,kappa,param_vec});
GradF = matlabFunction(GradS,'vars',{input_s,output_s,input_ns,output_ns,kappa,param_vec});
HessF = matlabFunction(HessS,'vars',{input_s,output_s,input_ns,output_ns,kappa,param_vec});
Pas1F = matlabFunction(Pas1S,'vars',{param_vec});
Pas2F = matlabFunction(Pas2S,'vars',{param_vec});
GrP1F = matlabFunction(GrP1S,'vars',{param_vec});
GrP2F = matlabFunction(GrP2S,'vars',{param_vec});

end

