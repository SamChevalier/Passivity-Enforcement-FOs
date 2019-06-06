function [Y_inf,fval] = PI_InferGen(data,basis_mats,enf_psvty,x0_opts,kappa)
% PI_INFERGEN:   Solve an admittance matrix inference problem.
%                Include no prior information, but regularize.
%
% Inputs:
% 1) data         Input data structure, containing two elements:
%                     1) data.ins     => 2x1 complex inputs (source frequnecy)
%                     2) data.outs    => 2x1 complex outputs (source frequnecy)
%                     3) data.inns    => 2x1 complex inputs (non-source frequnecy)
%                     4) data.outns   => 2x1 complex outputs (non-source frequnecy)
% 2) basis_mats   List of basis matrices which may be used in constructing
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
% 3) enf_psvty    This flag corresponds to passivity enforcement if
%                 enf_psvty = 1; otherwise, do not enforce passivity
% 4) x0_opts      Structure with x0 and trials
%                     1) x0_opts.x0     => starting position, if applicable
%                     2) x0_opts.TolFun => stopping tolerance
% 5) kappa        Regularizing constant
% 
% Outputs:
% 1) Y_inf        Inferred admittance matrix
% 2) fval         Associated objective function value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of included basis matrices
nb = sum(cell2mat(struct2cell(basis_mats)));

% Build cost, gradient and Hessian functions
[CostF,GradF,HessF,Pas1F,Pas2F,GrP1F,GrP2F] = PI_OptToolsGen(basis_mats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first three of these functions have the same five desired inputs:
%         1) input_s                      (2x1 complex)
%         2) output_s                     (2x1 complex)
%         3) input_ns                     (mx1 complex)
%         4) output_n                     (mx1 complex)
%         5) kappa regularizer            (1x1 real)
%         6) coefficient parameter values (nx1 real)
%                   ... where n is the number of assumed basis matrices,
%                       and "coefficient parameter values" are the
%                       variables which are being optimized over
%
% The following four only take the coefficient parameter values as inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymize
Cost = @(x)CostF(data.ins,data.outs,data.inns,data.outns,kappa,x);
Grad = @(x)GradF(data.ins,data.outs,data.inns,data.outns,kappa,x);
Hess = @(x)HessF(data.ins,data.outs,data.inns,data.outns,kappa,x);
Pas1 = @(x)Pas1F(x);
Pas2 = @(x)Pas2F(x);
GrP1 = @(x)GrP1F(x);
GrP2 = @(x)GrP2F(x);

% Set gradient and Hessian
options.SpecifyObjectiveGradient = true;
options.HessianFcn               = @hess_fun;

% Set stopping tolerance
options.OptimalityTolerance      = x0_opts.TolFun;

% Set passivity constraint if necessary
if enf_psvty == 1
    options.SpecifyConstraintGradient = true;
    const_fun = @ineq_con;
else
    const_fun = [];
end

% Starting point for initial run
if isempty(x0_opts.x0)
    x0 = zeros(nb,1);
else
    x0 = x0_opts.x0;
end

% Optimize
[Y_inf,fval,~,~] = fmincon(@opt_fun,x0,[],[],[],[],[],[],const_fun,options);

% Objective and gradient functions
function [f,g] = opt_fun(x)
    f = Cost(x);
    g = Grad(x);
end

% Hessian Function
function [h] = hess_fun(x)
    h = Hess(x);
end

% Passivity Constraint + Gradient
function [c,ceq,gc,gceq] = ineq_con(x)
    c    = [Pas1(x); Pas2(x)];
    ceq  = [];
    gc   = [GrP1(x); GrP2(x)];
    gceq = [];
end

end

