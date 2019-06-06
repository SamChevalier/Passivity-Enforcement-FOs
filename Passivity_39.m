%% 39 Bus System Test
clear variables; clc;

% Loads are given OU noise in the frequnecy domain, the generator inference
% problem is expanded, and current injections are defined with the same
% convention as shunt currents.

% ~ All generators are replaced with 3rd order + first order AVR ~

% Global Variables
global RSc DSc

% Damping and Resistance Factors
RSc = 1; % Nominal Resistance
DSc = 1; % Nominal Damping

%% EE's Rules
j  = sqrt(-1);
Md = [0 j; -j 0]; % Energy Function

% Frequnecy vector
f_vec = 1.95:0.01:2.05;
nf    = length(f_vec);

% Basis matrices
T1 = [1  0; 0  1]; jT1 = j*T1;
T2 = [1  0; 0 -1]; jT2 = j*T2;
T3 = [0  1; 1  0]; jT3 = j*T3;
T4 = [0 -1; 1  0]; jT4 = j*T4;

%% Initialize System
initpsat;                           % Initialize PSAT global variables
datafile        = 'NE39';           % Test case data file
runpsat(datafile,'data');           % Initialize datafile
Settings.freq   = 60;               % Change System Freq from default to 60
clpsat.readfile = 1;                % Read data file before running power flow
runpsat('pf');                      % Run power flow

% Parse Variables
Alg_Vars   = DAE.y;
State_Vars = DAE.x;

%% Parse Results
Va         = Alg_Vars(Bus.a);
Vm         = Alg_Vars(Bus.v);
Del        = State_Vars(Syn.delta);  % These will be replaced
Ef         = Alg_Vars(Syn.vf);       % These will be replaced

% Power Injection Results
Pgs = Bus.Pg;
Qgs = Bus.Qg;
Pls = Bus.Pl;
Qls = Bus.Ql;
Pi  = Pgs-Pls;
Qi  = Qgs-Qls;

% Determine Negative Current Injections
Ic = -conj((Pi+j*Qi)./(Vm.*exp(j*Va)));

%% Define new generator parameters (taken from Pai)
load('Gen_Params');
load('Load_Params');

% Random Generator Parameters
% % H    = 6*  (rand(10,1) + 0.5);
% % D    = 1*  (rand(10,1) + 0.5);
% % Xd   = 1.2*(rand(10,1) + 0.5);
% % Xdp  = 0.3*(rand(10,1) + 0.5);
% % Xq   = 1.1*(rand(10,1) + 0.5);
% % Td0p = 5*  (rand(10,1) + 0.5);
% % KA   = 20* (rand(10,1) + 0.5);
% % TA   = 0.5*(rand(10,1) + 0.5);

% Load Paramets
% % alpha_p = 2*rand(39,1);
% % alpha_q = 2*rand(39,1);
% % beta_p  = 2*rand(39,1);
% % beta_q  = 2*rand(39,1);

% Convert from (suspected) 13.8kV machine base to (suspected) 69kV base
Xd  = (13.8/69)^2*Xd;
Xdp = (13.8/69)^2*Xdp;
Xq  = (13.8/69)^2*Xq;

% Define parameter errors
H_e    = H.*(1+0.25*randn(10,1));
D_e    = D.*(1+0.25*randn(10,1));
Xd_e   = Xd.*(1+0.25*randn(10,1));
Xdp_e  = Xdp.*(1+0.25*randn(10,1));
Xq_e   = Xq.*(1+0.25*randn(10,1));
Td0p_e = Td0p.*(1+0.25*randn(10,1));
KA_e   = KA.*(1+0.25*randn(10,1));
TA_e   = TA.*(1+0.25*randn(10,1));

%% Solve for generator states based on new parameters
syms    d eqp Ef Vref
vars = [d eqp Ef Vref];

% Generator powers and voltages
Pgi = Pgs(30:39);
Qgi = Qgs(30:39);
V   = Vm(30:39);
T   = Va(30:39);

% Initialize Empty Vectors
del_vec  = zeros(10,1);
eqp_vec  = zeros(10,1);
Ef_vec   = zeros(10,1);
Vref_vec = zeros(10,1);

% Loop over generators
for ii = 1:10
    
    % Voltage and Current
    ed = V(ii)*sin(d-T(ii));
    eq = V(ii)*cos(d-T(ii));
    id = (eqp-eq)/Xdp(ii);
    iq = ed/Xq(ii);
    
    % Injections and dynamics
    eqs = [0       == (Ef - (Xd(ii)-Xdp(ii))*id-eqp)/Td0p(ii);
           0       == (KA(ii)*(Vref-V(ii))-Ef)/TA(ii);
           Pgi(ii) == ed*id + eq*iq;
           Qgi(ii) == id*eq - iq*ed];
    
    % Solution
    soln = vpasolve(eqs, vars, [NaN NaN; 0 4; 0 4; 0 4]);
    
    % Parse
    del_vec(ii)  = double(soln.d(1));
    eqp_vec(ii)  = double(soln.eqp(1));
    Ef_vec(ii)   = double(soln.Ef(1));
    Vref_vec(ii) = double(soln.Vref(1));
end

%% Injection
Inj_Bus  = 31;
Inj_Freq = 2;

% Representative Injection
Ir_inj   =  0.4547 + 0.0765i;
Ii_inj   = -1.7793 + 4.0604i;

% For a random injection
% % Ir_inj   = randn*exp(j*randn);
% % Ii_inj   = randn*exp(j*randn);

% Injection Magnitude
InjM = norm([Ir_inj; Ii_inj]);

% Injection Vector (Negative because original convection was opposite)
Js   = -[Ir_inj; Ii_inj];

%% Load Spectrum
tau = 1/(2*pi);
sig = 1e-3;

% Build FRF for filtering Gaussian noise
omega_loads = 2*pi*(f_vec);
FRF_OU      = 1./(j*omega_loads*tau + 1);
inj_ind     = find(f_vec == Inj_Freq);

%% Line Admittance Matrix (Diagonal)
nL   = Line.n;
Lf_i = Line.fr;
Lt_i = Line.to;
LR   = 2.5*Line.con(:,8);   % Line Resistance (increased)
LX   = Line.con(:,9);       % Line Reactance

% Build YL: Matrix of Line Admittance Matrices
YL     = zeros(2*nL,2*nL);
for ii = 1:nL
    YL((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = (1/(LR(ii)^2+LX(ii)^2))*[LR(ii)  LX(ii);
                                                                  -LX(ii)  LR(ii)];
end

%% Incidence Matrix
E  = zeros(2*nL,2*Bus.n);
Es = zeros(2*nL,2*Bus.n);
Er = zeros(2*nL,2*Bus.n);

% Loop over Lines
for ii = 1:nL
    line_from = Lf_i(ii);
    line_to   = Lt_i(ii);
    
    E((2*ii-1):(2*ii),(2*line_from-1):(2*line_from)) = eye(2,2);
    E((2*ii-1):(2*ii),(2*line_to-1):(2*line_to))     = -eye(2,2);
    
    % Sending End Voltage
    Es((2*ii-1):(2*ii),(2*line_from-1):(2*line_from)) = eye(2,2);
    
    % Receiving End Voltage
    Er((2*ii-1):(2*ii),(2*line_to-1):(2*line_to))     = eye(2,2);
end

%% Define Shunts: PQ Admittances
Vc = Vm.*exp(j*Va);
Vr = real(Vc);
Vi = imag(Vc);

% Power Injections
Pv_i         = zeros(Bus.n,1);
Pv_i(PQ.bus) = PQ.P0;         % Active Power Injection
Qv_i         = zeros(Bus.n,1);
Qv_i(PQ.bus) = PQ.Q0;         % Reactive Power Injection

% Use to compute the associated current injections
I_inj = conj((Pv_i+j*Qv_i)./Vc);
Ir    = real(I_inj);
Ii    = imag(I_inj);

% Power perturbations
S = zeros(39*2,nf);

% Load Buses
PQ_buses = PQ.con(:,1);

% Vector of currents injected by loads
I_inj_loads = zeros(2*39,nf);
V_RL        = zeros(2*39,nf);

% Matrix of PQ Admittances
Ym_PQ = zeros(2*Bus.n,2*Bus.n);
for ii = 1:Bus.n
    % Skip buses with no loads
    if ismember(ii,PQ_buses)
        % Skip the source bus
        if ii ~= Inj_Bus
            % Define Ya
            Ya = [alpha_p(ii)*Pv_i(ii)    beta_p(ii)*Pv_i(ii);
                  alpha_q(ii)*Qv_i(ii)    beta_q(ii)*Qv_i(ii)];
              
            % Define Ij
            Ij = [1 0; 0 j];
            
            % Define T1
            theta = Va(ii);
            vmag  = Vm(ii);
            T1m = [cos(theta) -vmag*sin(theta);
                   sin(theta)  vmag*cos(theta)];
               
             % Define AV & AI
             AV = [Vr(ii)   Vi(ii);
                   Vi(ii)  -Vr(ii)];
             AI = [Ir(ii)   Ii(ii);
                   -Ii(ii)  Ir(ii)];
               
             % Build Ypq
             ypq = AV\(Ya*Ij/T1m - AI);
             
             % Assemble
             Ym_PQ((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = ypq;
             
             % Also, define the current injection associated with the
             % loads' power perturbations
             for jj = 1:nf
                 u                   = FRF_OU(jj)*sig*exp(j*rand(1,1));
                 S((2*ii-1):2*ii,jj) = u.*[Pv_i(ii); Qv_i(ii)];
                     
                 % Transform to V (rectangular)
                 TM   = Ya*Ij/T1m;
                 
                 % If TM is rank 0, the current flow is 0
                 if rank(TM) < 2
                     V_RL((2*ii-1):2*ii,jj) = zeros(2,1);
                 else
                     V_RL((2*ii-1):2*ii,jj) = TM\S((2*ii-1):2*ii,jj);
                 end
                 
                 % Transform to I (rectangular)
                 I_RL = ypq*V_RL((2*ii-1):2*ii,jj);
                 
                 % Add to current injection vector - these currents flow
                 % *into* the loads
                 I_inj_loads((2*ii-1):2*ii,jj) = I_RL;
                 
% %                  % Test: if jj == inj_ind, no injection
% %                  if jj == inj_ind
% %                     I_inj_loads((2*ii-1):2*ii,jj) = 0*I_RL;
% %                  end
             end
        end
    end
end

%% Construct J Vectors
J = I_inj_loads;

% Add source injection
J((2*Inj_Bus-1):2*Inj_Bus,inj_ind) = J((2*Inj_Bus-1):2*Inj_Bus,inj_ind) + Js;

%% Define Shunts: Generator Admittances
Ym_Gen = zeros(2*Bus.n,2*Bus.n,nf);

% Call FRF function
[YaF] = YaF_3rd_AVR;

% Track eigenvalues and basis norms
eig_gens  = zeros(10,2);
bas_norms = zeros(10,1);

% Loop over relevant frequnecies
for jj = 1:nf
    % Omegas
    Omega_loop = omega_loads(jj);
    
    % Loop and build FRFs
    for ii  = 1:10
        kk  = ii + 29;
        
        % Skip the source bus (at forcing frequency
        if (kk == Inj_Bus) && (jj == inj_ind)
            % Arbitrary
        else
            FRF = YaF(D(ii),H(ii),KA(ii),Omega_loop,TA(ii),...
                Td0p(ii),Vi(ii),Vr(ii),Xd(ii),Xdp(ii), ...
                Xq(ii),del_vec(ii),eqp_vec(ii));
            
            % Assign to Shunt Matrix
            Ym_Gen((2*kk-1):(2*kk),(2*kk-1):(2*kk),jj) = FRF;
            
            % Track Eigenvalues of Matrices
            evs            = eig(FRF*Md + (FRF*Md)');
            eig_gens(ii,:) = evs(:);
            
            % Track basis matrix contributions
            [BasisCoeffs] = PI_BasisMats(FRF);
            bas_norms(ii) = norm(BasisCoeffs.aT1 + j*BasisCoeffs.bT1)/norm(FRF);
        end
    end
end

%% Generator Admittances: Errors for Initialization
Ym_Gen_Init = zeros(2*Bus.n,2*Bus.n);

% Loop and build FRFs
for ii  = 1:10
    kk  = ii + 29;
    
    % Don't skip the source bus
    FRF = YaF(D_e(ii),H_e(ii),KA_e(ii),2*pi*Inj_Freq,TA_e(ii),...
              Td0p_e(ii),Vi(ii),Vr(ii),Xd_e(ii),Xdp_e(ii), ...
              Xq_e(ii),del_vec(ii),eqp_vec(ii));
    
    % Assign to Shunt Matrix
    Ym_Gen_Init((2*kk-1):(2*kk),(2*kk-1):(2*kk)) = FRF;
end

%% Assemble "Ybus" Matrix at Each Frequnecy and Compute System Voltages
Vb        = zeros(39*2,nf);
Il        = zeros(46*2,nf);
Is_Gen    = zeros(39*2,1);
Is_Load   = zeros(39*2,1);
Is_Source = zeros(2,1);

% Loop
for jj = 1:nf
    
    % Call local generator matrix
    Ym_Gen_local = Ym_Gen(:,:,jj);
    Yb           = E'*YL*E + (Ym_Gen_local + Ym_PQ);
    
    % Solve the System
    Vb(:,jj) = -Yb\J(:,jj);
    
    % Shunt Currents?
    Is = (Ym_Gen_local + Ym_PQ)*Vb(:,jj);
    
    % Generator and Load Shunts
    Is_Gen(:,jj)         = (Ym_Gen_local)*Vb(:,jj);
    Is_Load(:,jj)        = (Ym_PQ)*Vb(:,jj);
    Is_Source(:,inj_ind) = Js;
    
    % Line Currents
    Il(:,jj) = YL*E*Vb(:,jj);
    
    % Build DEF Matrix
    M  = zeros(2*Bus.n,2*Bus.n);
    for ii = 1:Bus.n
        M((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Md;
    end
    
    % Mulitplication Matrix - Buses
    K     = zeros(Bus.n,2*Bus.n);
    kk    = 1;
    for ii = 1:Bus.n
        K(ii,kk)   = 1;
        K(ii,kk+1) = 1;
        kk = kk+2;
    end
    
    % Shunt power
    D_source = K*real(conj(diag(Vb(:,jj)))*M*J(:,jj));
    
    % Quantify DEF - matches if Resistance = 0
    D_shunts = K*real(conj(diag(Vb(:,jj)))*M*Is);
end

%% Compute Actual Power Perturbations at Loads
S_act   = zeros(39*2,1);
Ild_act = I_inj_loads + Is_Load;

for ii = 1:Bus.n
    % Skip buses with no loads
    if ismember(ii,PQ_buses)
        % Skip the source bus
        if ii ~= Inj_Bus
            
            % Define T1
            theta = Va(ii);
            vmag  = Vm(ii);
            T1m = [cos(theta) -vmag*sin(theta);
                   sin(theta)  vmag*cos(theta)];
            
            % Define AV & AI
            AV = [Vr(ii)   Vi(ii);
                  Vi(ii)  -Vr(ii)];
            AI = [Ir(ii)   Ii(ii);
                  -Ii(ii)  Ir(ii)];
            
            % Define Actual Power
            S_act( (2*ii-1):2*ii ) = AI*Vb((2*ii-1):2*ii,inj_ind) + AV*Ild_act((2*ii-1):2*ii,inj_ind);
        end
    end
end

%% Quantify DEF on Lines - Forcing Frequnecy
M  = zeros(2*nL,2*nL);
for ii = 1:nL
    M((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Md;
end

% Quantify
D_vecS = conj(diag(Il(:,inj_ind)))*M*Es*Vb(:,inj_ind);
D_vecR = conj(diag(Il(:,inj_ind)))*M*Er*Vb(:,inj_ind);

% Mulitplication Matrix - Lines
K  = zeros(Line.n,2*Line.n);
kk = 1;
for ii = 1:Line.n
    K(ii,kk)   = 1;
    K(ii,kk+1) = 1;
    kk = kk+2;
end

% Final DEF
Df1S = real(K*D_vecS);
Df1R = real(K*D_vecR);

% Remove source generator from the shunt matrix
Ym_GenNS                                                      = Ym_Gen(:,:,inj_ind);
Ym_GenNS((2*Inj_Bus-1):(2*Inj_Bus),(2*Inj_Bus-1):(2*Inj_Bus)) = zeros(2,2);

% Build the Dynamic Ward Equivalent 
Yb     = E'*YL*E + (Ym_GenNS + Ym_PQ);
[Yeq]  = PI_Dyn_Ward_Equiv(Yb,Inj_Bus);
Vs     = Vb((2*Inj_Bus-1):(2*Inj_Bus),inj_ind);

% Test source energy from Y
D_sourceY = sum(real(Vs'*Md*(Yeq*Vs)));

% Test DWE eigenvalues
[V,Deig] = eig(Yeq*Md + (Yeq*Md)');

%% Find Injection
% % Yeqi   = inv(Yeq);
% % for ii = 1:10000
% %     Jsl           = randn(2,1) + j*randn(2,1);
% %     D_sourceL(ii) = sum(real((Yeqi*Jsl)'*Md*(Jsl)));
% %     Jsl_L(:,ii)   = Jsl;
% % end

%% Draw System Map
clf
[C] = Draw_39_Map(Df1S,Df1R);
% export_fig Full_Sys.png -m9

%% Solve Load Inference - Just at Forcing Frequnecy

% Solve at loads
Ym_PQinf = zeros(Bus.n*2,Bus.n*2);
YaS      = zeros(2,2,39);

for ii = [1:29 39];
    % Define inputs/outputs
    Vrec_pert = Vb( (2*ii-1):(2*ii) , inj_ind );
    
    % Transform to polar (define T1m)
    theta = Va(ii);
    vmag  = Vm(ii);
    T1m = [cos(theta) -vmag*sin(theta);
           sin(theta)  vmag*cos(theta)];
    Vpolar_pert = T1m\Vrec_pert;
    Ij          = [1 0; 0 j];
    
    % Data
    y = S_act((2*ii-1):(2*ii));  % Loads are ZIP
    u = Ij*Vpolar_pert;          % Input
    
    data.in  = u;
    data.out = y;
    
    % Define use-able basis matrices
    basis_mats.T1  = 1;
    basis_mats.T2  = 1;
    basis_mats.T3  = 1;
    basis_mats.T4  = 1;
    basis_mats.jT1 = 0;
    basis_mats.jT2 = 0;
    basis_mats.jT3 = 0;
    basis_mats.jT4 = 0;
    
    % Passivity
    enf_psvty = 0;
    
    % Initial conditions?
    x0_opts.x0     = [];
    x0_opts.trials = 1;
    x0_opts.TolFun = 1e-20;

    % Prior model?
    prior.Y0       = [];
    prior.lambda   = 0;

    % Solve inference problem
    [Y_inf,Y_infs,fvals] = PI_Infer(data,basis_mats,enf_psvty,x0_opts,prior);
    
    % Use Y_inf result to build admittance
    YaS(:,:,ii) = Y_inf(1)*T1 + Y_inf(2)*T2 + Y_inf(3)*T3 + Y_inf(4)*T4;
    
    % Define AV & AI
    AV = [Vr(ii)   Vi(ii);
          Vi(ii)  -Vr(ii)];
    AI = [Ir(ii)   Ii(ii);
         -Ii(ii)   Ir(ii)];

    % Populate PQ inference matrix
    Ym_PQinf((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = AV\(YaS(:,:,ii)*Ij/T1m - AI);
end

%% Solve Generator Inference
Ym_Geninf          = zeros(Bus.n*2,Bus.n*2);
n_inj_ind          = 1:nf;
n_inj_ind(inj_ind) = [];

% Loop over *all* gens
for ii = 30:39;
    
    % Input/Output
    us = Vb((2*ii-1):2*ii,inj_ind);          % Input (forcing frequnecy)
    ys = Is_Gen((2*ii-1):2*ii,inj_ind);      % Output (forcing frequnecy)
    
    % Input/Output
    uns = Vb((2*ii-1):2*ii,n_inj_ind);     % Input (non-forcing frequnecies)
    yns = Is_Gen((2*ii-1):2*ii,n_inj_ind); % Output (non-forcing frequnecies)
    
    % Reshape to column
    uns = reshape(uns,[2*nf-2 1]); % Doesn't include source frequnecy
    yns = reshape(yns,[2*nf-2 1]);
    
    % Use injection if source bus
    if ii == Inj_Bus
        ys = Js;
    end
    
    % Data
    data.ins   = us;
    data.outs  = ys;
    data.inns  = uns;
    data.outns = yns;
    
    % Define use-able basis matrices
    basis_mats.T1  = 0;
    basis_mats.T2  = 1;
    basis_mats.T3  = 1;
    basis_mats.T4  = 1;
    basis_mats.jT1 = 0;
    basis_mats.jT2 = 1;
    basis_mats.jT3 = 1;
    basis_mats.jT4 = 1;
    
    % Passivity
    enf_psvty = 1;
    
    % Initial conditions? Set to random
    ygt               = Ym_Gen_Init((2*ii-1):(2*ii),(2*ii-1):(2*ii));
    [BasisCoeffs]     = PI_BasisMats(ygt);
    x0_opts.x0        = cell2mat(struct2cell(BasisCoeffs));
    x0_opts.x0([1 5]) = [];
    x0_opts.x0        = [];
    x0_opts.TolFun    = 1e-20;
    
    % Regularizer
    kappa = 1e6;

    % Solve
    [Y_inf,fvals] = PI_InferGen(data,basis_mats,enf_psvty,x0_opts,kappa);
    
    % Build Yg
    Yg = Y_inf(1)*T2  +   Y_inf(2)*T3 +   Y_inf(3)*T4 + ...
       j*Y_inf(4)*T2  + j*Y_inf(5)*T3 + j*Y_inf(6)*T4;
     
    % Populate Gen inference matrix
    Ym_Geninf((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Yg;
end

%% Perform Kron Reductions at Gens - Original Simple Testing
error_Y = zeros(10,1);    % Admittance matrix error
DWE_vec = zeros(2,2,10);
error   = zeros(10,1);    % Simple current prediction error
kk      = 1;

% Measured generator currents
Is_Gen_M                            = Is_Gen(:,inj_ind);
Is_Gen_M((2*Inj_Bus-1):(2*Inj_Bus)) = Js;

% Loop over all generators
for ii = 30:39
    % Remove suspected source gen
    Ym_Geninf_RmvSrc                                  = Ym_Geninf;
    Ym_Geninf_RmvSrc((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = zeros(2,2);
    
    % Remove suspected source load - not necessary actually
    Ym_PQinf_RmvSrc = Ym_PQinf;
    
    % Build network model
    Yb = E'*YL*E + (Ym_Geninf_RmvSrc + Ym_PQinf_RmvSrc);
    
    % Build DWE
    bus             = ii;
    DWE_vec(:,:,kk) = -PI_Dyn_Ward_Equiv(Yb,bus);
    
    % Compute error associated with DWE estimate - normalized?
    error(ii-29) = norm(Is_Gen_M((2*ii-1):(2*ii))-DWE_vec(:,:,kk)*Vb((2*ii-1):(2*ii),inj_ind))/norm(Is_Gen_M((2*ii-1):(2*ii)));
    
    % Compute error via admittance differentials
    error_Y(ii-29) = norm(DWE_vec(:,:,kk) - Ym_Geninf((2*ii-1):(2*ii),(2*ii-1):(2*ii)))/norm(DWE_vec(:,:,kk));
    
    % Increment
    kk = kk+1;
    
end

%% % New Test
%  Loop over all generators, and loop again, removing one generator each
%  time. Record the error for each removal
ll = 1;
aa = 1;
bb = 1;
cc = 1;
dd = 1;
ee = 1;
ff = 1;
gg = 1;
hh = 1;
pp = 1;
qq = 1;
error_sub = zeros(9,10);

% Loop over all generators
for ii = 30:39
    bus_vec = 30:39;
    
    % Remove ii
    bus_vec(bus_vec == ii) = [];
    
    % Loop over remaining
    kk = 1;
    for jj = bus_vec
        inj_bus  = ii;
        edge_bus = jj;
        
        % Remove source and edge generators
        Ym_Geninf_RmvSrc = Ym_Geninf;
        Ym_Geninf_RmvSrc((2*inj_bus-1):(2*inj_bus),(2*inj_bus-1):(2*inj_bus))     = zeros(2,2);
        Ym_Geninf_RmvSrc((2*edge_bus-1):(2*edge_bus),(2*edge_bus-1):(2*edge_bus)) = zeros(2,2);
        
        % Remove source and edge generators - no need to remove source load
        Ym_PQinf_RmvSrc = Ym_PQinf;
        %Ym_PQinf_RmvSrc((2*inj_bus-1):(2*inj_bus),(2*inj_bus-1):(2*inj_bus))     = zeros(2,2);
        %Ym_PQinf_RmvSrc((2*edge_bus-1):(2*edge_bus),(2*edge_bus-1):(2*edge_bus)) = zeros(2,2);
        
        % Build network model
        Yb = E'*YL*E + (Ym_Geninf_RmvSrc + Ym_PQinf_RmvSrc);
        
        % Determine submatrices
        [m1,m2,m3,m4,~,~] = PI_Ward_Sub_Mats(Yb,inj_bus,edge_bus);
        
        % Compute
        Vs = Vb((2*inj_bus-1):(2*inj_bus),inj_ind);
        Ie = Is_Gen_M((2*edge_bus-1):(2*edge_bus));
        Is = Is_Gen_M((2*inj_bus-1):(2*inj_bus));
        Ip = (m1 - m2*(m4\m3))*Vs - m2*(m4\Ie);
        
        % Save the effects that removing generator "jj" has
        if jj == 30
            er30(aa) = norm(Is+Ip);
            aa=aa+1;
        elseif jj == 31
            er31(bb) = norm(Is+Ip);
            bb=bb+1;
        elseif jj == 32
            er32(cc) = norm(Is+Ip);
            cc=cc+1;
        elseif jj == 33
            er33(dd) = norm(Is+Ip);
            dd=dd+1;
        elseif jj == 34
            er34(ee) = norm(Is+Ip);
            ee=ee+1;
        elseif jj == 35
            er35(ff) = norm(Is+Ip);
            ff=ff+1;
        elseif jj == 36
            er36(gg) = norm(Is+Ip);
            gg=gg+1;
        elseif jj == 37
            er37(hh) = norm(Is+Ip);
            hh=hh+1;
        elseif jj == 38
            er38(pp) = norm(Is+Ip);
            pp=pp+1;
        elseif jj == 39
            er39(qq) = norm(Is+Ip);
            qq = qq+1;
        end
            
        % Full error
        error_sub(kk,ll) = norm(Is+Ip);
        kk = kk+1;
    end
    ll = ll+1;
end

%% Test Removal Method
inj_bus  = 31;
edge_bus = 39;

% Weave vectors
vinj                          = (31*2-1):(2*31);
vedge                         = (35*2-1):(2*35);
Ym_Geninf_RmvSrc              = Ym_Geninf;
Ym_Geninf_RmvSrc(vinj,vinj)   = zeros(2,2);
Ym_Geninf_RmvSrc(vedge,vedge) = zeros(2,2);

% Remove suspected source load - not necessary
Ym_PQinf_RmvSrc              = Ym_PQ;
% % Ym_PQinf_RmvSrc(vinj,vinj)   = zeros(2,2);
% % Ym_PQinf_RmvSrc(vedge,vedge) = zeros(2,2);

% Build network model
Yb                        = E'*YL*E + (Ym_Geninf_RmvSrc + Ym_PQinf_RmvSrc);
[m1,m2,m3,m4,Ybn,bus_vec] = PI_Ward_Sub_Mats(Yb,inj_bus,edge_bus);

% Compute
Vs = Vb(vinj,inj_ind);
Ie = Is_Gen_M((2*edge_bus-1):(2*edge_bus));
Is = Is_Gen_M(vinj);
Ip = (m1 - m2*(m4\m3))*Vs - m2*(m4\Ie);


%% Plot Error - Original
clf
gen_inds = [30 32:39];
semilogy(gen_inds,error([1 3:end]),'xblack','LineWidth',1.1,'MarkerSize',7)
hold on
semilogy(31,error(2),'*red','LineWidth',1.1,'MarkerSize',7)
hold on
for ii = 30:39
    line([ii ii],[error(ii-29) 10e-6],'color','black')
end
xlim([29 40])
ylim([1e-4 1e2])
set(gca,'XTick',30:39)
set(gca,'XTicklabel',30:39)
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1])
set(gca,'YTicklabel',{'','10^{-3}','','10^{-1}','','10^{1}'})
set(gca,'FontName','Times','FontSize',14)
ylabel('$\rm{Prediction \;\, Error}$','Interpreter','latex','FontSize',14)
xlabel('$\rm{Generator \;\, Bus}$','Interpreter','latex','FontSize',14)
legend({'${\rm Non-Source \; Gen}$','${\rm Source \; Gen}$'},'Interpreter','latex','box','off','Location','NW','FontSize',13)
set(gca,'YGrid','on');
set(gca,'YMinorGrid','on')
set(gcf,'Units','inches','Position',[0 0 9 2.5])
tightfig(gcf); % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% Plot Generator Removals
c30 = [0, 0.4470, 0.7410];
c31 = [0.8500, 0.3250, 0.0980];
c32 = [0.9290, 0.6940, 0.1250];
c33 = [0.4940, 0.1840, 0.5560];
c34 = [0.4660, 0.6740, 0.1880];
c35 = [0.3010, 0.7450, 0.9330];
c36 = [0.6350, 0.0780, 0.1840];
c37 = [0, 0.75, 0.75];
c38 = [0, 0.5, 0];
c39 = [0 0 0];

% Build Semilogy Plot
clf
semilogy(35)
hold on
a1  = scatter(1,1,100,c30,'LineWidth',1.5,'Marker','*');
a2  = scatter(1,1,100,c31,'LineWidth',1.5,'Marker','h');
a3  = scatter(1,1,100,c32,'LineWidth',1.5,'Marker','*');
a4  = scatter(1,1,100,c33,'LineWidth',1.5,'Marker','*');
a5  = scatter(1,1,100,c34,'LineWidth',1.5,'Marker','*');
a6  = scatter(1,1,100,c35,'LineWidth',1.5,'Marker','*');
a7  = scatter(1,1,100,c36,'LineWidth',1.5,'Marker','*');
a8  = scatter(1,1,100,c37,'LineWidth',1.5,'Marker','*');
a9  = scatter(1,1,100,c38,'LineWidth',1.5,'Marker','*');
a10 = scatter(1,1,100,c39,'LineWidth',1.5,'Marker','*');

% Add legend now
legend([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10],{'${\rm Gen \; 30}$','${\rm Gen \; 31}$','${\rm Gen \; 32}$','${\rm Gen \; 33}$','${\rm Gen \; 34}$','${\rm Gen \; 35}$','${\rm Gen \; 36}$','${\rm Gen \; 37}$','${\rm Gen \; 38}$','${\rm Gen \; 39}$',},'Interpreter','latex','box','off','Location','SE','FontSize',13)
text(39.43,1.3,'${\rm Gen \; Eliminated\!\!:}$','Interpreter','latex','FontSize',13)

% Plot Data
scatter(30*ones(8,1),error_sub(2:end,1),100,[c32; c33; c34; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(31*ones(9,1),error_sub(:,2),100,[c30; c32; c33; c34; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(32*ones(8,1),error_sub([1 3:end],3),100,[c30; c33; c34; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(33*ones(8,1),error_sub([1 3:end],4),100,[c30; c32; c34; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(34*ones(8,1),error_sub([1 3:end],5),100,[c30; c32; c33; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(35*ones(8,1),error_sub([1 3:end],6),100,[c30; c32; c33; c34; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(36*ones(8,1),error_sub([1 3:end],7),100,[c30; c32; c33; c34; c35; c37; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(37*ones(8,1),error_sub([1 3:end],8),100,[c30; c32; c33; c34; c35; c36; c38; c39],'LineWidth',1.5,'Marker','*')
scatter(38*ones(8,1),error_sub([1 3:end],9),100,[c30; c32; c33; c34; c35; c36; c37; c39],'LineWidth',1.5,'Marker','*')
scatter(39*ones(8,1),error_sub([1 3:end],10),100,[c30; c32; c33; c34; c35; c36; c37; c38],'LineWidth',1.5,'Marker','*')

% Plot Generator 31 Data with star markers
scatter(30,error_sub(1,1),100,c31,'LineWidth',2,'Marker','h')
scatter(32,error_sub(2,3),100,c31,'LineWidth',2,'Marker','h')
scatter(33,error_sub(2,4),100,c31,'LineWidth',2,'Marker','h')
scatter(34,error_sub(2,5),100,c31,'LineWidth',2,'Marker','h')
scatter(35,error_sub(2,6),100,c31,'LineWidth',2,'Marker','h')
scatter(36,error_sub(2,7),100,c31,'LineWidth',2,'Marker','h')
scatter(37,error_sub(2,8),100,c31,'LineWidth',2,'Marker','h')
scatter(38,error_sub(2,9),100,c31,'LineWidth',2,'Marker','h')
scatter(39,error_sub(2,10),100,c31,'LineWidth',2,'Marker','h')

% Connect scatter plot points with lines
xlim([29.5 42])
ylim([4e-6 0.4e1])
plot([31:39],er30,'color',c30,'LineWidth',1)
plot([30 32:39],er31,'color',c31,'LineWidth',1)
plot([30:31 33:39],er32,'color',c32,'LineWidth',1)
plot([30:32 34:39],er33,'color',c33,'LineWidth',1)
plot([30:33 35:39],er34,'color',c34,'LineWidth',1)
plot([30:34 36:39],er35,'color',c35,'LineWidth',1)
plot([30:35 37:39],er36,'color',c36,'LineWidth',1)
plot([30:36 38:39],er37,'color',c37,'LineWidth',1)
plot([30:37 39],er38,'color',c38,'LineWidth',1)
plot([30:38],er39,'color',c39,'LineWidth',1)

% Highlight 2 stories
color_val = [0 0 0];
face_color_val = [c30 0.1];
rectangle('Position',[29.6 16^(-4) 9.8 0.00017],'LineWidth' ,1,'EdgeColor',color_val,'Curvature',[0.2,0.5],'LineStyle','-','facecolor',face_color_val);
face_color_val = [c31 0.1];
rectangle('Position',[30.5 45^(-3) 1 0.00065],'LineWidth' ,1,'EdgeColor',color_val,'Curvature',[0.3,0.2],'LineStyle','-','facecolor',face_color_val);

% Add labels etc.
xlabel({'${\rm Generator \; Index}$'},'Interpreter','latex','FontSize',14)
ylabel({'${\rm Current \; Prediction \; Error}$'},'Interpreter','latex','FontSize',14)
text(29.6,16.2^-4,'${\rm a}$','Interpreter','latex','FontSize',13)
text(30.25,13^-3,'${\rm b}$','Interpreter','latex','FontSize',13)
set(gca,'XTick',[30 31 32 33 34 35 36 37 38 39])
set(gca,'XTicklabel',{'30', '31', '32', '33', '34', '35', '36', '37', '38', '39'})
set(gca,'FontName','Times','FontSize',14)
set(gcf,'Units','inches','Position',[0 0 9 3.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"