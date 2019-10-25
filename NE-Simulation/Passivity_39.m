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
[V,Deig] = eig(Yeq*Md + (Yeq*Md)')
stop = 1;
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