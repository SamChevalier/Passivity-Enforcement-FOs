function [Y_inf,Y_infs,fvals] = PI_Infer(data,basis_mats,enf_psvty,x0_opts,prior)
% PI_INFER:       Solve an admittance matrix inference problem. Interpret
%                 the "prior_mats" inclusions/exclusions as constraints.
%
% Inputs:
% 1) data         Input data structure, containing two elements:
%                     1) data.in    => 2x1 complex inputs
%                     2) data.out   => 2x1 complex outputs
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
%                     2) x0_opts.trials => number of random trials, if applicable
%                     3) x0_opts.TolFun => stopping tolerance
% 5) prior        Structure with prior data
%                     1) prior.Y0       => numerical prior matrix
%                     2) prior.lambda   => numerical lambda (regularizer)
% 
% Outputs:
% 1) Y_inf        Inferred admittance matrix (average if many)
% 2) Y_infs       3D matrix of all solutions
% 3) fvals        Associated objective function values
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of included basis matrices
nb = sum(cell2mat(struct2cell(basis_mats)));

% Deconstruct the prior matrix (if supplied)
if isempty(prior.Y0)
    prior_coeffs = zeros(nb,1);
    prior.lambda = 0;      % Define just in case
else
    % Remove basis matrices not included in the "basis_mats" structure
    [prior_coeffs] = PI_BasisMats(prior.Y0);
    prior_coeffs   = cell2mat(struct2cell(prior_coeffs));
    prior_coeffs(~cell2mat(struct2cell(basis_mats))) = [];
end

% Build cost, gradient and Hessian functions
[CostF,GradF,HessF,Pas1F,Pas2F,GrP1F,GrP2F] = PI_OptTools(basis_mats);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first three of these functions have the same five desired inputs:
%         1) input                        (2x1 complex)
%         2) output                       (2x1 complex)
%         3) lambda regularizer           (1x1 real)
%         4) prior vector                 (nx1 real)
%         5) coefficient parameter values (nx1 real)
%                   ... where n is the number of assumed basis matrices,
%                       and "coefficient parameter values" are the
%                       variables which are being optimized over
%
% The following four only take the coefficient parameter values as inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anonymize
Cost = @(x)CostF(data.in,data.out,prior.lambda,prior_coeffs,x);
Grad = @(x)GradF(data.in,data.out,prior.lambda,prior_coeffs,x);
Hess = @(x)HessF(data.in,data.out,prior.lambda,prior_coeffs,x);
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
[Y_inf,fval,exitflag,output] = fmincon(@opt_fun,x0,[],[],[],[],[],[],const_fun,options);

% Iterations?
if isempty(x0_opts.trials) || (x0_opts.trials==1)
    Y_infs = Y_inf;
    fvals  = fval;
else
    % Loop and optimize
    Y_infs = zeros(nb,x0_opts.trials);
    fvals  = zeros(x0_opts.trials,1);
    
    for ii = 1:x0_opts.trials
        x0           = 10*randn(nb,1);
        [Y_inf,fval] = fmincon(@opt_fun,x0,[],[],[],[],[],[],const_fun,options);
        Y_infs(:,ii) = Y_inf;
        fvals(ii)    = fval;
    end
    
    % Average
    Y_inf = mean(Y_infs,2);
end

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

