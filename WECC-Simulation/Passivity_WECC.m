%% Use PSAT to Simulate Case 1 (Single Source = Bus 4 @ 0.86Hz)
%  OTF Source: http://web.eecs.utk.edu/~kaisun/Oscillation/1f.html
%  Sam Chevalier
clear variables; clc;

% Resolved issue: 1) If a voltage magnitude of the PF solution is
% outside the specified voltage range, loads become impedances, 2)
% don't specify impedance characteristics in the perturbation file, and 3)
% don't use COI! It will fail you, master.

% EE's rule
j = sqrt(-1);

% Source variables
global Source_V Amplitude_V Frequency_V vref0 Fl

% Load params
global alpha_p alpha_q beta_p beta_q ij

% PSAT load/integration globals
global intstep VP_nos Vl0v VQ_nos tcorr P0 Q0 PQ_std Varout Vlb

%%% Careful! All bus indices are shifted by 1. Therfore, "Bus 4" %%%
%%% (source Bus) is actually Bus #3 when calling its variables.  %%%

%% 1) * * * Basis Matrices * * *
Tb1 = [1  0; 0  1]; jTb1 = j*Tb1;
Tb2 = [1  0; 0 -1]; jTb2 = j*Tb2;
Tb3 = [0  1; 1  0]; jTb3 = j*Tb3;
Tb4 = [0 -1; 1  0]; jTb4 = j*Tb4;

% Initialize first time
alpha_p = 0; alpha_q = 0; beta_p = 0; beta_q = 0;
WECC_179

%% 1A) Initial Data

% Set Load Parameters
nLoads  = length(PQ.con);
alpha_p = 2*ones(nLoads,1); % Constant R Load
alpha_q = 2*ones(nLoads,1); % Constant X Load
beta_p  = 1*ones(nLoads,1);
beta_q  = 1*ones(nLoads,1);

% Intialize system
WECC_179
npq = size(PQ.con,1);           % Number of load buses

% Load Indices
load_inds = PQ.con(:,1)-1;

%% 1B) Solve Power Flow
initpsat;
datafile = 'WECC_179';
runpsat(datafile,'data');
Settings.freq   = 60;
Settings.maxvar = 5000;    % Increase # of Variables
runpsat('pf');             % Run the Almighty Power Flow %

P0  = Fl.con(:,2);   % Initial value of loads' active power
Q0  = Fl.con(:,5);   % Initial value of loads' reactive power
Tf  = Fl.con(1,8);   % Filter time constant

% Parse Power Flow Results Data
Angles   = DAE.y(Bus.a);
Voltages = DAE.y(Bus.v);
EMFs     = Syn.vf0;

% Determine the base voltages of loads
Vl0v = Voltages(load_inds);

% Store 3rd Order Generator Initialization Data
i1        = Syn.e1q;
Syn3.eqp  = DAE.x(i1(1));
Syn3.Ef   = DAE.x(Exc.vf);

%% 1C) Oscillation Source Data
%  Sources are generator indices

% Exciter Reference
Source_Bus  = 4-1; % Index, not bus number
Source_V    = 1;
Amplitude_V = 0.10;
Frequency_V = 0.86 + .1;
vref0       = Exc.vref0(Source_V);

% Band around forcing frequency to use for inference
epF = 0.03;

%% 1D) Initialize Time Domain Simulation
intstep = 0.0333333;
tcorr   = 1;
PQ_std  = 1e-3;
tbegin  = 0;
tfinal  = 120;
runpsat('Vref_perturb','pert');     % Perturbation file
Settings.freq   = 60;               % Change System Freq from default to 60
clpsat.readfile = 1;                % Read data file before running power flow
VP_nos  = zeros(npq,1e5);           % Vector of noise for load P
VQ_nos  = zeros(npq,1e5);           % Vector of noise for load Q
ij      = 1;                        % Initialize noise index

% SETTINGS FOR TIME DOMAIN SIMULATION
Settings.coi   = 0;                % Do *NOT* use center of inertia for synchronous machines
Settings.t0    = tbegin;           % Initial simulation time
Settings.tf    = tfinal;           % Final simulation time
Settings.pq2z  = 0;                % Do not convert PQ loads to constant impedance (Z) loads after power flow
Settings.fixt  = 1;                % Enable fixed time-step solver
Settings.tstep = intstep;          % Fixed time step value
nL             = Line.n + Ltc.n + Phs.n + Hvdc.n + Lines.n; % Num circuit elements
Varname.idx    = 1:DAE.n + DAE.m + 2*Bus.n + 6*nL + 2;      % Plot Variable Indexes (ask for current to be included)

% Bus Voltage Magnitude Indices
ix_Vm  = DAE.n+Bus.n+1:DAE.n+2*Bus.n;    
Vlb    = ix_Vm(load_inds);

% Run Time Domain Simulation
runpsat('td');

%% 1E) Define Output Variables
ix_Va  = DAE.n+1:DAE.n+Bus.n;                                              % Index of voltage angles
ix_Vm  = DAE.n+Bus.n+1:DAE.n+2*Bus.n;                                      % Index of voltage magnitudes
ix_Iij = DAE.n + DAE.m + 2*Bus.n + 4*nL+1:DAE.n + DAE.m + 2*Bus.n + 5*nL;  % Index of line currents
ix_Iji = DAE.n + DAE.m + 2*Bus.n + 5*nL+1:DAE.n + DAE.m + 2*Bus.n + 6*nL;  % Index of line currents
ix_Pij = DAE.n + DAE.m + 2*Bus.n + 1:DAE.n + DAE.m + 2*Bus.n + nL;         % Index of line active power flows (i -> j)
ix_Qij = DAE.n + DAE.m + 2*Bus.n +2*nL +1:DAE.n + DAE.m + 2*Bus.n + 3*nL;  % Index of line reactive powers flows (i -> j)
ix_de  = Syn.delta;                                                        % Index of rotor angles
ix_om  = Syn.omega;                                                        % Index of generator speed
ix_Pi  = DAE.n + DAE.m + 1:DAE.n + DAE.m + Bus.n;                          % Index of Active Power Injections
ix_Qi  = DAE.n + DAE.m + Bus.n + 1:DAE.n + DAE.m + 2*Bus.n;                % Index of Reactive Power Injections
ix_Pm  = DAE.n + Syn.pm;  

% Call all of the algebraic variables saved during the TD Simulation
state_ind  = 1:DAE.n;
alg_ind    = (DAE.n + 1):(DAE.n + DAE.m);
Alg_Vars   = Varout.vars(:,alg_ind);
State_Vars = Varout.vars(:,state_ind);

% Output variable values
Va        = Varout.vars(:,ix_Va);   % Va: Bus voltage angles
Vm        = Varout.vars(:,ix_Vm);   % Vm: Bus voltage magnitudes
time      = Varout.t;
deltas    = Varout.vars(:,ix_de);
omegas    = Varout.vars(:,ix_om);
Il_ij     = Varout.vars(:,ix_Iij);
Il_ji     = Varout.vars(:,ix_Iji);
Pij       = Varout.vars(:,ix_Pij);  % Active Power Flows
Qij       = Varout.vars(:,ix_Qij);  % Reactive Power Flows
P         = Varout.vars(:,ix_Pi);   % Active Power Injections
Q         = Varout.vars(:,ix_Qi);   % Reactive Power Injections
Pm        = Varout.vars(:,ix_Pm);   % Mechanical Powers
dw        = Alg_Vars(:,Fl.dw);

% Other Data
V_ref   = Alg_Vars(:,Exc.vref);
eqp_ind = Syn.e1q;
eqp     = State_Vars(:,eqp_ind(1));

% Trim Noise Vectors
VP_nos = VP_nos(:,length(time));
VQ_nos = VQ_nos(:,length(time));

%% 1F) Save all Voltage and Current Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save('Sim_Data','Vm','Va','deltas','time','V_ref','eqp','EMFs','intstep','Syn3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1G) What are the Active Power Injections at Each Generator Terminal?
Gen_inds = Syn.con(:,1);
P_Gens   = zeros(size(P,1),length(Gen_inds));
for ii = 1:length(Gen_inds)
    Gen_bus      = Gen_inds(ii) - 1; % Because Bus #4 is the 3rd bus
    P_Gens(:,ii) = detrend(P(:,Gen_bus),'constant');
end

%% %%%%%%%%%%%%%%%%%%%%% Part 2: Test Results %%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2A) Build Complex Currents, Voltages, and Power
nb = Bus.n;
Vc = zeros(size(Vm));
Sc = zeros(size(Vm));
Ic = zeros(size(Vm));

for ii = 1:nb
    Vc(:,ii) = Vm(:,ii).*exp(j*unwrap(Va(:,ii)));
    Sc(:,ii) = P(:,ii) + j*Q(:,ii);
    Ic(:,ii) = conj(-Sc(:,ii)./Vc(:,ii)); %- j*Vc(:,ii)*Line.con();
end

%% 2B) Trim First 20s of PMU Data & (FFT) Process
nds  = round(20/intstep);
r.Vc = Vc(nds:end,:);
r.Sc = Sc(nds:end,:);
r.Ic = Ic(nds:end,:);

% Initialize TD
r.Vr = zeros(size(r.Vc));
r.Vi = zeros(size(r.Vc));
r.Ir = zeros(size(r.Vc));
r.Ii = zeros(size(r.Vc));

% Initialize FD
n      = length(r.Vr);
[f,~]  = Apply_FFT_N(r.Vr(:,1),intstep,n);
r.y_Vr = zeros(length(f),nb);
r.y_Vi = zeros(length(f),nb);
r.y_Ir = zeros(length(f),nb);
r.y_Ii = zeros(length(f),nb);

% Decompose
for ii = 1:nb
    r.Vr(:,ii) = real(r.Vc(:,ii));
    r.Vi(:,ii) = imag(r.Vc(:,ii));
    r.Ir(:,ii) = real(r.Ic(:,ii));
    r.Ii(:,ii) = imag(r.Ic(:,ii));
    
    % Apply FFT
    n          = length(r.Vr);
    HWnd       = hann(n);
    [f,y_Vr] = Apply_FFT_N(HWnd.*detrend(r.Vr(:,ii),'constant'),intstep,n);
    [~,y_Vi] = Apply_FFT_N(HWnd.*detrend(r.Vi(:,ii),'constant'),intstep,n);
    [~,y_Ir] = Apply_FFT_N(HWnd.*detrend(r.Ir(:,ii),'constant'),intstep,n);
    [~,y_Ii] = Apply_FFT_N(HWnd.*detrend(r.Ii(:,ii),'constant'),intstep,n);
    
    % FFT vectors
    r.f          = f;
    r.y_Vr(:,ii) = y_Vr;
    r.y_Vi(:,ii) = y_Vi;
    r.y_Ir(:,ii) = y_Ir;
    r.y_Ii(:,ii) = y_Ii;
end

% Build Full Vector
Vb = zeros(2*nb,length(f));
Ib = zeros(2*nb,length(f));
for kk = 1:length(f)
    for ii = 1:nb
        Vb(2*ii-1:2*ii,kk) = [r.y_Vr(kk,ii); r.y_Vi(kk,ii)];
        Ib(2*ii-1:2*ii,kk) = [r.y_Ir(kk,ii); r.y_Ii(kk,ii)];
    end
end

% Index closest to FO frequency
[~,ind] = min(abs(f - Frequency_V));

%% 2C) Build Line Admittance Matrix (Diagonal)
nL   = Line.n;
Lf_i = Line.fr;
Lt_i = Line.to;
LR   = Line.con(:,8);   % Line Resistance
LX   = Line.con(:,9);   % Line Reactance
Ybus = Line.Y;          % PSAT's Ybus   

% Build YL: Matrix of Line Admittance Matrices
v      = zeros(nL,2);
kk     = 1;
YL     = zeros(2*nL,2*nL);
for ii = 1:nL
    lf = Lf_i(ii);
    lt = Lt_i(ii);
    
    % Frist, test if we have encountered this line already
    vecF = lf == v(:,1);
    vecT = lt == v(:,2);
    
    if sum(2 == vecF + vecT)
        % Skip this line - it is a duplicate
    else
        y  = Ybus(lf,lt);
        g  = -real(y);
        b  = imag(y);
        YL((2*kk-1):(2*kk),(2*kk-1):(2*kk)) = [g  b
                                              -b  g];
        % The value of "b" is negative in the bottom left because we use the
        % convention 1/(R+jX) = G+jB
        
        % Add to some vector
        v(kk,:) = [lf lt];
        kk      = kk+1;
    end
end

% Number of distinct lines
nL_dist   = kk-1;
Line_list = v;
YL        = YL(1:2*nL_dist,1:2*nL_dist);

%% 2D) Build Incidence Matrix
E  = zeros(2*nL_dist,2*Bus.n);
Es = zeros(2*nL_dist,2*Bus.n);
Er = zeros(2*nL_dist,2*Bus.n);

% Loop over Lines, but skip repeats
for ii = 1:nL_dist
    line_from = Line_list(ii,1);
    line_to   = Line_list(ii,2);
    
    E((2*ii-1):(2*ii),(2*line_from-1):(2*line_from)) = eye(2,2);
    E((2*ii-1):(2*ii),(2*line_to-1):(2*line_to))     = -eye(2,2);
    
    % Sending End Voltage
    Es((2*ii-1):(2*ii),(2*line_from-1):(2*line_from)) = eye(2,2);
    
    % Receiving End Voltage
    Er((2*ii-1):(2*ii),(2*line_to-1):(2*line_to))     = eye(2,2);
end

%% 2E) Define ZIP Load Shunts
Ym_PQ  = zeros(2*Bus.n,2*Bus.n);
Ym_PQp = zeros(2*Bus.n,2*Bus.n); % Polar power representation

% Load Buses (mind the "-1" applied to the index)
PQ_buses = Fl.con(:,1)-1;

jj = 1;
for ii = PQ_buses'
    
    % Define T1
    theta = Angles(ii);
    vmag  = Voltages(ii);
    T1m = [cos(theta) -vmag*sin(theta);
           sin(theta)  vmag*cos(theta)];
    
    % Compute SS currents/voltages associated with loads
    IrSS = real(conj( (P0(jj)+j*Q0(jj))/(vmag*exp(j*theta)) ));
    IiSS = imag(conj( (P0(jj)+j*Q0(jj))/(vmag*exp(j*theta)) ));
    VrSS = real(vmag*exp(j*theta));
    ViSS = imag(vmag*exp(j*theta));
    
    Ya = [(1/vmag)*alpha_p(jj)*P0(jj)    beta_p(jj)*P0(jj)/(2*pi*60/Tf);
          (1/vmag)*alpha_q(jj)*Q0(jj)    beta_q(jj)*P0(jj)/(2*pi*60/Tf)];
    
    % Define Ij
    Ij = [1 0; 0 j];
    
    % Define AV & AI
    AV = [VrSS   ViSS;
          ViSS  -VrSS];
    AI = [IrSS   IiSS;
         -IiSS   IrSS];
    
    % Build Ypq
    ypq = AV\(Ya*Ij/T1m - AI);
    
    % Assemble
    Ym_PQ((2*ii-1):(2*ii),(2*ii-1):(2*ii))  = ypq;
    Ym_PQp((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Ya;
    
    % Increment
    jj = jj+1;
end 

%% 2F) Build Shunts - Line and Non-Line
Ym_Shunt = zeros(2*Bus.n,2*Bus.n);
Ysh      = sparse(Shunt.con(:,1)-1,Shunt.con(:,1)-1,Shunt.con(:,5)+j*Shunt.con(:,6),nb,nb);
Ybus_S   = full(Ybus + Ysh);

% Shunt Buses
for ii = 1:nb
    ys = sum(Ybus_S(ii,:));
    gs = real(ys);
    bs = imag(ys);
    Ys = [gs -bs;
          bs  gs];

    % Populate Shunt Matrix
    Ym_Shunt((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Ys;
end

%% 2Ga) Build Generator Shunts
Ym_Gen = zeros(2*Bus.n,2*Bus.n);

% 2nd Order Machine Parameter Vectors
GIv    = Syn.con(2:end,1);
GIvD   = Syn.con(2:end,1)-1;
Mv     = Syn.con(2:end,18);
Dv     = Syn.con(2:end,19);
Xdpv   = Syn.con(2:end,9);
EMFv   = EMFs(2:end);
n2nd   = length(EMFv);
Vmv    = Vm(:,GIvD);
Vav    = Va(:,GIvD);
Deltav = deltas(:,2:end);

% 3rd Order Generator (Source 1)
G3     = Syn.con(1,1);
GD3    = G3 - 1;
Vm3    = Vm(:,GD3);
Va3    = Va(:,GD3);
Delta3 = deltas(:,1);
M3     = Syn.con(1,18);
D3     = Syn.con(1,19);
xd3    = Syn.con(1,8);
xdp3   = Syn.con(1,9);
xq3    = Syn.con(1,13);
Td0p3  = Syn.con(1,11);

% Exciter Data
KA3    = Exc.con(1,5);
TA3    = Exc.con(1,6);

% Initialized Values: Neither of these affect the Admittance Matrix
Vr3 = V_ref(1);
Pm3 = 0;

%% 2Gb) 2nd Order Generator
syms del_a w_a V_a T_a Ef_a D_Ya M_Ya X_Ya Pm_a Omega
X   = [del_a w_a];
Uv  = [V_a T_a];
w0    = 2*pi*60;
Pe    = V_a*Ef_a*sin(del_a-T_a)/X_Ya;
F_vec = [w0*w_a;
        (Pm_a - Pe - D_Ya*w_a)/M_Ya];
G_vec = [(1/X_Ya)*sqrt(V_a^2 + Ef_a^2 - 2*Ef_a*V_a*cos(T_a-del_a)); 
         atan((V_a*sin(T_a)-Ef_a*sin(del_a))/(V_a*cos(T_a)-Ef_a*cos(del_a)))-pi];
JacFX  = jacobian(F_vec,X);
JacFUv = jacobian(F_vec,Uv);
JacGX  = jacobian(G_vec,X);
JacGUv = jacobian(G_vec,Uv);
Imgen  = zeros(n2nd,1);
Phgen  = zeros(n2nd,1);
Y2nd_cell = cell(n2nd,1);
for ii = 1:n2nd    % Loop Over 2nd Order Gens
    Ef_a     = EMFv(ii);
    D_Ya     = Dv(ii);
    M_Ya     = Mv(ii);
    X_Ya     = Xdpv(ii);
    del_a    = Deltav(1,ii);
    w_a      = 0;
    V_a      = Vmv(1,ii);
    T_a      = Vav(1,ii);
    JacFX_n  = subs(JacFX);
    JacFUv_n = subs(JacFUv);
    JacGX_n  = subs(JacGX);
    JacGUv_n = subs(JacGUv);
    Y_2nd    = (JacGX_n/(j*Omega*eye(length(X)) - JacFX_n))*JacFUv_n + JacGUv_n;
    Y2nd_cell{ii} = matlabFunction(Y_2nd);
    
    % Collect Current Data
    I         = (V_a*exp(j*T_a) - Ef_a*exp(j*del_a))/(j*X_Ya);
    Imgen(ii) = abs(I);
    Phgen(ii) = angle(I);
end

%% 2Gc) 3rd Order Gen
syms del_a w_a eqp_a Ef_a V_a T_a Vr_a Omega
X   = [del_a w_a eqp_a Ef_a];
Uv  = [V_a T_a];
Ue  = [Vr_a];
ed = V_a*sin(del_a-T_a);
eq = V_a*cos(del_a-T_a);
id = (eqp_a-eq)/xdp3;
iq = (ed)/xq3;
Pe = ed*id + eq*iq;
F_vec = [w0*w_a;
        (Pm3 - Pe - w_a)/(M3);
        (Ef_a - (xd3-xdp3)*id - eqp_a)/Td0p3;
        (KA3*(Vr3-V_a)-Ef_a)/(TA3)];
Ir    = -sin(del_a)*id - cos(del_a)*iq;
Ii    =  cos(del_a)*id - sin(del_a)*iq;
G_vec = [sqrt(id^2+iq^2); angle(Ir+j*Ii)];
JacFX  = jacobian(F_vec,X);
JacFUv = jacobian(F_vec,Uv);
% % JacFUe = jacobian(F_vec,Ue);
JacGX  = jacobian(G_vec,X);
JacGUv = jacobian(G_vec,Uv);
del_a  = Delta3(1);
w_a    = 0;
eqp_a  = Syn3.eqp;
Ef_a   = Syn3.Ef;
V_a    = Vm3(1);
T_a    = Va3(1);
Vr_a   = Vr3(1);
JacFX_n  = subs(JacFX);
JacFUv_n = subs(JacFUv);
JacGX_n  = subs(JacGX);
JacGUv_n = subs(JacGUv);
Y3rd     = matlabFunction((JacGX_n/(j*Omega*eye(length(X)) - JacFX_n))*JacFUv_n + JacGUv_n);

%% 2Gd) Populate FRF Matrix
Gen2_buses = Syn.con(2:end,1)-1;
GenS       = Syn.con(1,1)-1;
jj         = 1;
for ii = Gen2_buses'
    Omega_n = 2*pi*(Frequency_V);
    Y_gen   = Y2nd_cell{jj}(Omega_n);
    
    % The admittance must be converted to cartesian
    theta = Angles(ii);
    vmag  = Voltages(ii);
    T1v   = [cos(theta) -vmag*sin(theta);
             sin(theta)  vmag*cos(theta)];
    Imag = Imgen(jj);
    phi  = Phgen(jj);
    T1i  = [cos(phi) -Imag*sin(phi);
            sin(phi)  Imag*cos(phi)];
    
    Ym_Gen((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = T1i*Y_gen/T1v;
    
    % Increment
    jj = jj+1;
end

% 3rd Order FRF (Apply Transformation)
Y3    = Y3rd(Omega_n);
theta = Angles(Source_Bus);
vmag  = Voltages(Source_Bus);
Tv    = [cos(theta) -vmag*sin(theta);
         sin(theta)  vmag*cos(theta)];
PgS   = DAE.y(Syn.p(Source_V)); % Call P
QgS   = DAE.y(Syn.q(Source_V)); % Call Q
Ig3   = conj((PgS+j*QgS)/(vmag*exp(j*theta)));
Imag  = abs(Ig3);
phi   = angle(Ig3);
Ti    = [cos(phi) -Imag*sin(phi);
         sin(phi)  Imag*cos(phi)];
Ym_Gen((2*GenS-1):(2*GenS),(2*GenS-1):(2*GenS)) = Ti*Y3/Tv;

%% 2H) Assemble Full Yb Matrix
Yb_full  = E'*YL*E + (Ym_Gen + Ym_PQ + Ym_Shunt);

% Assemble Yb Matrix at Source Terminals
Yb_Sh = (Ym_Gen + Ym_PQ + Ym_Shunt);
Yb_Sh(Source_Bus*2-1:Source_Bus*2,Source_Bus*2-1:Source_Bus*2) = zeros(2,2);
YbS   = E'*YL*E + Yb_Sh;

% Reduce System and Check Eigenvalues
[Yeq]  = PI_Dyn_Ward_Equiv(YbS,Source_Bus);
M      = [0 j; -j 0];
0.5*eig(M*Yeq + (M*Yeq)')

% Shunts Injections
Ish = (E'*YL*E)*Vb(:,ind);

% Compute the injections at all nodes
jj  = 1;
Def = zeros(nb,1);
for ii = 1:nb
    Def(jj) = real(Vb(ii*2-1:ii*2,ind)'*M*Ish(ii*2-1:ii*2));
    jj      = jj+1;
end

% Test accuracy of results
Test_Inj  = (E'*YL*E + Ym_Gen + Ym_PQ + Ym_Shunt)*Vb(:,ind);

% Test again with the following :)
% Test_Inj2 = (E'*YL*E + Ym_Gen + Ym_PQ + Ym_Shunt)*Vb_eR;

%% %%%%%%%%%%%%%%%%%%%%% Part 3: Solve Inference %%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3A) Get Load Powers
% Load Buses (mind the "-1" applied to the index)
PQ_buses = Fl.con(:,1)-1;
rPQ.P    = zeros(length(Vm),nLoads);
rPQ.Q    = zeros(length(Vm),nLoads);

jj = 1;
for ii = PQ_buses'
    % Positive mean flowing *into* the bus
    rPQ.P(:,jj) = P0(jj).*(1+VP_nos(jj,:)').*(Vm(:,ii)./Vm(1,ii)).^alpha_p(jj).*(1+dw(:,jj)).^beta_p(jj);
    rPQ.Q(:,jj) = Q0(jj).*(1+VQ_nos(jj,:)').*(Vm(:,ii)./Vm(1,ii)).^alpha_q(jj).*(1+dw(:,jj)).^beta_q(jj);
    jj = jj+1;
end

% Trim the vectors
rPQT.P = rPQ.P(nds:end,:);
rPQT.Q = rPQ.Q(nds:end,:);

%% 3B) Get Generator Currents
id = (eqp - Vm3.*cos(Delta3-Va3) )/xdp3;
iq = (Vm3.*sin(Delta3-Va3))/xq3;

% Convert 3rd Order Gen Currents to Polar
r3.Im  = zeros(size(id,1),1);
r3.Ia  = zeros(size(id,1),1);
for ii = 1:length(id)
    T = -[ sin(Delta3(ii))  cos(Delta3(ii));
          -cos(Delta3(ii))  sin(Delta3(ii))];
    I = T*[id(ii); iq(ii)];
    r3.Im(ii) = abs(I(1) + j*I(2));
    r3.Ia(ii) = angle(I(1) + j*I(2));
end

% Unwrap 4th Order Generator
r3.Ia = unwrap(r3.Ia);

% 2nd Order Gens
r2.Im  = zeros(length(Vm),n2nd);
r2.Ia  = zeros(length(Vm),n2nd);
jj     = 1;
for ii = GIvD'    % Loop Over 2nd Order Gens
    ig = (Vc(:,ii) - EMFv(jj)*exp(j*Deltav(:,jj)))/(j*Xdpv(jj));
    r2.Im(:,jj) = abs(ig);
    r2.Ia(:,jj) = unwrap(angle(ig));
    jj = jj+1;
end

% Trim the vectors
r.Vm   = Vm(nds:end,:);
r.Va   = Va(nds:end,:);
r3T.Im = r3.Im(nds:end);
r3T.Ia = r3.Ia(nds:end);
r2T.Im = r2.Im(nds:end,:);
r2T.Ia = r2.Ia(nds:end,:);

%% 3C) Add Measurement Noise

%% 3D) Take FFT of all Inference Data

% Loop over voltage
[f,~]  = Apply_FFT_N(r.Vm(:,1),intstep,length(r.Vm(:,1)));
r.y_Vm = zeros(length(f),Bus.n);
r.y_Va = zeros(length(f),Bus.n);
for ii = 1:Bus.n
    n            = length(r.Vm(:,ii));
    HWnd         = hann(n);
    [f,y_Vm]     = Apply_FFT_N(HWnd.*detrend(r.Vm(:,ii),'constant'),intstep,n);
    [~,y_Va]     = Apply_FFT_N(HWnd.*detrend(r.Va(:,ii),'constant'),intstep,n);
    r.y_Vm(:,ii) = y_Vm;
    r.y_Va(:,ii) = y_Va;
end

% Loop over PQ loads
r.y_P = zeros(length(f),nLoads);
r.y_Q = zeros(length(f),nLoads);
for ii = 1:nLoads
    n           = length(r.Vm(:,1));
    HWnd        = hann(n);
    [f,y_P]     = Apply_FFT_N(HWnd.*detrend(rPQT.P(:,ii),'constant'),intstep,n);
    [~,y_Q]     = Apply_FFT_N(HWnd.*detrend(rPQT.Q(:,ii),'constant'),intstep,n);
    r.y_P(:,ii) = y_P;
    r.y_Q(:,ii) = y_Q;
end

% Loop over 2nd order generator currents
r.y_Im2 = zeros(length(f),n2nd);
r.y_Ia2 = zeros(length(f),n2nd);
for ii = 1:n2nd
    n           = length(r.Vm(:,1));
    HWnd        = hann(n);
    [f,y_Im]    = Apply_FFT_N(HWnd.*detrend(r2T.Im(:,ii),'constant'),intstep,n);
    [~,y_Ia]    = Apply_FFT_N(HWnd.*detrend(r2T.Ia(:,ii),'constant'),intstep,n);
    r.y_Im2(:,ii) = y_Im;
    r.y_Ia2(:,ii) = y_Ia;
end

% 3rd order generator currents
[f,y_Im] = Apply_FFT_N(HWnd.*detrend(r3T.Im,'constant'),intstep,n);
[~,y_Ia] = Apply_FFT_N(HWnd.*detrend(r3T.Ia,'constant'),intstep,n);
r.y_Im3  = y_Im;
r.y_Ia3  = y_Ia;

% Set fft indices
[~,indFO]      = min(abs(f - Frequency_V));
[~,indFO_PepF] = min(abs(f - Frequency_V+epF));
[~,indFO_MepF] = min(abs(f - Frequency_V-epF));

% Set up vector of indices with no FO frequency
nFO_vec                   = indFO_PepF:indFO_MepF;
nFO_vec(nFO_vec == indFO) = [];

%% 3E) Solve Load Inference

% Solve at loads
Ym_PQinf = zeros(2*Bus.n,2*Bus.n);

% Define usable basis matrices
basis_mats.T1  = 1;
basis_mats.T2  = 1;
basis_mats.T3  = 1;
basis_mats.T4  = 1;
basis_mats.jT1 = 0;
basis_mats.jT2 = 0;
basis_mats.jT3 = 0;
basis_mats.jT4 = 0;

% Build cost, gradient and Hessian functions
[PQfuncs] = PI_OptTools_Fast(basis_mats);

% Loop over all loads
jj = 1;
for ii = PQ_buses'
    % Define inputs/outputs
    Vp = [r.y_Vm(indFO,ii);
          r.y_Va(indFO,ii)];
    Ij = [1 0; 0 j]; 
    u  = Ij*Vp;                % Input
    y  = [r.y_P(indFO,jj);     % Output
          r.y_Q(indFO,jj)];
    
    % Define data structure
    data.in  = u;
    data.out = y;
      
    % Initial conditions
    x0_opts.x0     = [];
    x0_opts.TolFun = 1e-14;
    
    % Prior model?
    prior.Y0       = [];
    prior.lambda   = 0;
    
    % Solve inference problem
    [Y_inf,~] = PI_Infer_Fast(data,basis_mats,x0_opts,prior,PQfuncs);
    
    % Use Y_inf result to build admittance
    Ya = Y_inf(1)*Tb1 + Y_inf(2)*Tb2 + Y_inf(3)*Tb3 + Y_inf(4)*Tb4;
    
    % Define T1
    theta = Angles(ii);
    vmag  = Voltages(ii);
    T1m   = [cos(theta) -vmag*sin(theta);
             sin(theta)  vmag*cos(theta)];
    
    % Compute SS currents/voltages associated with loads
    IrSS = real(conj( (P0(jj)+j*Q0(jj))/(vmag*exp(j*theta)) ));
    IiSS = imag(conj( (P0(jj)+j*Q0(jj))/(vmag*exp(j*theta)) ));
    VrSS = real(vmag*exp(j*theta));
    ViSS = imag(vmag*exp(j*theta));
    
    % Define AV & AI
    AV = [VrSS   ViSS;
          ViSS  -VrSS];
    AI = [IrSS   IiSS;
         -IiSS   IrSS];

    % Populate PQ inference matrix
    Ym_PQinf((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = AV\(Ya*Ij/T1m - AI);
    
    % Increment
    jj = jj+1;
    
end

%% 3F) Solve Generator Inference
Gen_buses = Syn.con(:,1)-1;

% Solve at loads
Ym_Geninf = zeros(2*Bus.n,2*Bus.n);

% Define usable basis matrices
basis_mats.T1  = 0;
basis_mats.T2  = 1;
basis_mats.T3  = 1;
basis_mats.T4  = 1;
basis_mats.jT1 = 0;
basis_mats.jT2 = 1;
basis_mats.jT3 = 1;
basis_mats.jT4 = 1;

% Build cost, gradient and Hessian functions
[Gfuncs] = PI_OptToolsGen_Fast(basis_mats);

% Loop over all generators
jj = 1;
for ii = Gen_buses(1:end)'
    
    % First, build the Transformation Matrices
    theta = Angles(ii);
    vmag  = Voltages(ii);
    T1v   = [cos(theta) -vmag*sin(theta);
             sin(theta)  vmag*cos(theta)];
    
    % Input & output at forcing frequnecy
    us = [r.y_Vm(indFO,ii);
          r.y_Va(indFO,ii)];
    
    % Input & output at NON-forcing frequnecy
    unsVM = r.y_Vm(nFO_vec,ii);
    unsVA = r.y_Va(nFO_vec,ii);
    uns   = zeros(2*length(unsVM),1);
    uns(1:2:end) = unsVM;
    uns(2:2:end) = unsVA;
    
    % Repeat for currents (differentiate between 2nd and 3rd order gens)
    if jj == 1
        ys = [r.y_Im3(indFO);
              r.y_Ia3(indFO)];
        ynsIM = r.y_Im3(nFO_vec);
        ynsIA = r.y_Ia3(nFO_vec);
        yns   = zeros(2*length(ynsIM),1);
        yns(1:2:end) = ynsIM;
        yns(2:2:end) = ynsIA;
        
        % Generator current transformation
        Imag = r3.Im(1);
        phi  = r3.Ia(1);
        T1i  = [cos(phi) -Imag*sin(phi);
                sin(phi)  Imag*cos(phi)];
    else
        ys = [r.y_Im2(indFO,jj-1);
              r.y_Ia2(indFO,jj-1)];
        ynsIM = r.y_Im2(nFO_vec,jj-1);
        ynsIA = r.y_Ia2(nFO_vec,jj-1);
        yns   = zeros(2*length(ynsIM),1);
        yns(1:2:end) = ynsIM;
        yns(2:2:end) = ynsIA;
        
        % Generator current transformation
        Imag = Imgen(jj-1); % Subtract 1 becasue we're looping over the
        phi  = Phgen(jj-1); % full set of generators, not just 2nd order
        T1i  = [cos(phi) -Imag*sin(phi);
                sin(phi)  Imag*cos(phi)];
    end
    
    % % % Transform ALL inputs and outputs % % %
    us = T1v*us;
    ys = T1i*ys;
    for kk = 1:length(nFO_vec)
        yns(2*kk-1:2*kk) = T1i*yns(2*kk-1:2*kk);
        uns(2*kk-1:2*kk) = T1v*uns(2*kk-1:2*kk);
    end
    
    % Data
    data.ins   = us;
    data.outs  = ys;
    data.inns  = uns;
    data.outns = yns;
    
    % Passivity
    enf_psvty = 1;
    
    % Initial conditions? Set to random
    x0_opts.x0     = [];
    x0_opts.TolFun = 1e-15;
    
    % Regularizer
    kappa = 1e2;

    % Solve
    [Y_inf,fvals] = PI_InferGen_Fast(data,basis_mats,enf_psvty,x0_opts,kappa,Gfuncs);
    
    % Build Yg
    Yg = Y_inf(1)*Tb2  +   Y_inf(2)*Tb3 +   Y_inf(3)*Tb4 + ...
       j*Y_inf(4)*Tb2  + j*Y_inf(5)*Tb3 + j*Y_inf(6)*Tb4;
     
    % Populate Gen inference matrix
    Ym_Geninf((2*ii-1):(2*ii),(2*ii-1):(2*ii)) = Yg;
    
    % Increment
    jj = jj+1;
end

%% 3G) Transform all Bus Voltages & Gen Currents to "Effective" Rectangular
Vb_eR  = zeros(2*Bus.n,1);
Ig_eR  = zeros(2*(n2nd+1),1);
Ig_eRs = zeros(2*Bus.n,1); % Full system

% Loop over bus voltages
for ii = 1:Bus.n
    theta = Angles(ii);
    vmag  = Voltages(ii);
    T1v   = [cos(theta) -vmag*sin(theta);
             sin(theta)  vmag*cos(theta)];
    Vb_eR(2*ii-1:2*ii) = T1v*[r.y_Vm(indFO,ii);
                              r.y_Va(indFO,ii)];
end

% Loop over generator currents
jj = 1;
for ii = Gen_buses'
    if jj == 1   % 3rd order
        % Generator current transformation
        Imag = r3.Im(1);
        phi  = r3.Ia(1);
        T1i  = [cos(phi) -Imag*sin(phi);
                sin(phi)  Imag*cos(phi)];
        Ig_eR(2*jj-1:2*jj) = T1i*[r.y_Im3(indFO);
                                  r.y_Ia3(indFO)];
                              
        % Also, put the current in a full system-wide vector
        Ig_eRs(2*ii-1:2*ii) = Ig_eR(2*jj-1:2*jj);
        
    else % 2nd order
        % Generator current transformation
        Imag = Imgen(jj-1); % Subtract 1 becasue we're looping over the
        phi  = Phgen(jj-1); % full set of generators, not just 2nd order
        T1i  = [cos(phi) -Imag*sin(phi);
                sin(phi)  Imag*cos(phi)];
        Ig_eR(2*jj-1:2*jj) = T1i*[r.y_Im2(indFO,jj-1);
                                  r.y_Ia2(indFO,jj-1)];
        
        % Also, put the current in a full system-wide vector
        Ig_eRs(2*ii-1:2*ii) = Ig_eR(2*jj-1:2*jj);
    end
    
    % Increment
    jj = jj+1;
end

% Test again with the following :)
Test_Inj2 = (E'*YL*E + Ym_Gen + Ym_PQ + Ym_Shunt)*Vb_eR;

% Shunts Injections
Ish_eR = (E'*YL*E)*Vb_eR;

% Compute the injections at all nodes
jj  = 1;
Def_eR = zeros(nb,1);
for ii = 1:nb
    Def_eR(jj) = real(Vb_eR(ii*2-1:ii*2)'*M*Ish_eR(ii*2-1:ii*2));
    jj         = jj+1;
end

%% 3H) Generator Removal Test
% Loop over all generators, and loop again, removing one generator each
% time. Record the error for each removal.
error_mat2 = zeros(200,200);
ind_vec    = ones(200,1);

% Loop over all generators
for ii = Gen_buses'
    bus_vec = Gen_buses;
    
    % Remove ii
    bus_vec(bus_vec == ii) = [];
    
    % Loop over remaining
    kk = 1;
    for jj = bus_vec'
        inj_bus  = ii;  % Suspected source
        edge_bus = jj;  % Generator to remove
        
        % Remove source and edge generators
        Ym_Geninf_RmvSrcs = Ym_Geninf;
        Ym_Geninf_RmvSrcs((2*inj_bus-1):(2*inj_bus),(2*inj_bus-1):(2*inj_bus))     = zeros(2,2);
        Ym_Geninf_RmvSrcs((2*edge_bus-1):(2*edge_bus),(2*edge_bus-1):(2*edge_bus)) = zeros(2,2);
        
        % Build network model
        Yb = (E'*YL*E + Ym_Shunt) + (Ym_Geninf_RmvSrcs + Ym_PQinf);
        
        % Determine submatrices
        [m1,m2,m3,m4,~,~] = PI_Ward_Sub_Mats(Yb,inj_bus,edge_bus);
        
        % Compute
        Vs = Vb_eR((2*inj_bus-1):(2*inj_bus));
        Ie = Ig_eRs((2*edge_bus-1):(2*edge_bus));
        Is = Ig_eRs((2*inj_bus-1):(2*inj_bus));
        Ip = (m1 - m2*(m4\m3))*Vs - m2*(m4\Ie);
        % There is a negative here because of the current convention
        
        % Save the effects that removing generator "jj" has
        % error1 = norm(Is+Ip);
        error2 = norm(Is+Ip)/norm(Is);
        % error3 = norm(Is+Ip)/(2*norm(Is-Ip));
        
        % Record error
        hh                = ind_vec(jj);  % call the index
        error_mat2(jj,hh) = error2;       % store the data
        ind_vec(jj)       = hh+1;         % increase the increment
        
    end
end

% figure
% title('e1')
% semilogy(error_mat1(Gen_buses,1:29)','*')

figure
title('e2')
semilogy(error_mat2(Gen_buses,1:29)','*')

% figure
% title('e3')
% semilogy(error_mat3(Gen_buses,1:29)','*')

%% Paper Plot 1) Plot DEF
clf
DEF_vars = abs(Def_eR(Gen_buses)); % one has value of -epsilon
semilogy(1,DEF_vars(1),'*','markersize',12,'color','red','LineWidth',1.5)
legend({'${\rm Source \; Gen}$'},'Interpreter','latex','box','off','Location','NE','FontSize',13)
hold on
semilogy(2:29,DEF_vars(2:end),'x','markersize',8,'color','black','LineWidth',1.25)
for ii = 1:length(Gen_buses)
    line([ii; ii],[1e-15; DEF_vars(ii)],'Color','black')
end
grid on
ylim([1e-8 1e-3])
set(gca,'FontName','Times','FontSize',14)
xlabel('${\rm Generator \; Index}$','Interpreter','latex','FontSize',14)
ylabel('${\rm Dissipating \; Power \;}P^{\star}$','Interpreter','latex','FontSize',14)

set(gcf,'Units','inches','Position',[0 0 9 2.5])
tightfig(gcf); % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% Paper Plot 2) Plot Generator Removals
load('error_mat2_Master')
nG = length(Syn.con);
sc = [0.8500, 0.3250, 0.0980];
Cv = [0      0.4470 0.7410;
      0.9290 0.6940 0.1250;
      0.4940 0.1840 0.5560;
      0.4660 0.6740 0.1880;
      0.3010 0.7450 0.9330;
      0.6350 0.0780 0.1840;
      0      0.75   0.75;
      0      0.5    0;
      0      0      0];

% Parse Master
close all
semilogy(35)
hold on

emM = error_mat2_Master(Gen_buses,1:28,1);
for ii = 1:nG
    bl = 1:nG;
    bl(ii) = [];
    if ii == 1
        % Plot last
    else
        jj = ii - 9*floor(ii/9);
        if jj == 0
            jj = ii/9;
        end
        scatter(bl,100*emM(ii,:),100,Cv(jj,:),'LineWidth',1.5,'Marker','*');
        plot(bl,100*emM(ii,:),'color',Cv(jj,:),'LineWidth',0.5)
    end
end

% Source
plot(2:nG,100*emM(1,:),'color',sc,'LineWidth',2.5)
scatter(2:nG,100*emM(1,:),150,sc,'LineWidth',2,'Marker','h');

% Plot characteristics
xlim([-0.5 29.5])
ylim([2e-1 1000])

color_val      = [0 0 0];
face_color_val = [sc 0.1];
rectangle('Position',[0.4 4*10^(-1) 1.15 7.5],'LineWidth' ,1,'EdgeColor',color_val,'Curvature',[0.3,0.2],'LineStyle','-','facecolor',face_color_val);

% Add legend now
a = scatter(35,50,100,sc,'LineWidth',1.1,'Marker','h');
b = scatter(35,50,100,'black','LineWidth',1.1,'Marker','*');
legend([a,b],{'${\rm Source \; Gen}$','${\rm Other \; Gens}$'},'Interpreter','latex','box','off','Location','SE','FontSize',12.5)
text(24,3,'${\rm Gen \; Eliminated\!\!:}$','Interpreter','latex','FontSize',12.5)

xlabel({'${\rm Generator \; Index}$'},'Interpreter','latex','FontSize',14)
ylabel({'${\rm Current \; Prediction \; Error \; }(\%)$'},'Interpreter','latex','FontSize',14)
text(3.9,2.6,'${\rm a}$','Interpreter','latex','FontSize',14)
text(0.1,10.1,'${\rm b}$','Interpreter','latex','FontSize',14)
set(gca,'XTick',1:29)
set(gca,'XTicklabel',{'1', '', '', '', '5', '', '', '', '', '10', '', '', '', '', '15', '', '', '', '', '20', '', '', '', '', '25', '', '', '', '29'})
set(gca,'YTicklabel',{'1', '10', '100', '1000'})
set(gca,'FontName','Times','FontSize',14)
set(gcf,'Units','inches','Position',[0 0 9 3.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%%
% 1) Plot Scatter of Source
scatter(30*ones(8,1),error_sub(2:end,1),100,[c32; c33; c34; c35; c36; c37; c38; c39],'LineWidth',1.5,'Marker','*')

a2  = scatter(1,1,100,c31,'LineWidth',1.5,'Marker','h');

% c30 = [0, 0.4470, 0.7410];
% c31 = [0.8500, 0.3250, 0.0980];
% c32 = [0.9290, 0.6940, 0.1250];
% c33 = [0.4940, 0.1840, 0.5560];
% c34 = [0.4660, 0.6740, 0.1880];
% c35 = [0.3010, 0.7450, 0.9330];
% c36 = [0.6350, 0.0780, 0.1840];
% c37 = [0, 0.75, 0.75];
% c38 = [0, 0.5, 0];
% c39 = [0 0 0];

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
legend([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10],{'${\rm Gen \; 30}$','${\rm Gen \; 31}$','${\rm Gen \; 32}$','${\rm Gen \; 33}$','${\rm Gen \; 34}$','${\rm Gen \; 35}$','${\rm Gen \; 36}$','${\rm Gen \; 37}$','${\rm Gen \; 38}$','${\rm Gen \; 39}$',},'Interpreter','latex','box','off','Location','SE','FontSize',12)
text(39.43,1.3,'${\rm Gen \; Eliminated\!\!:}$','Interpreter','latex','FontSize',12.5)

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

%% Plot Current Injection Vector J
Jv = abs(Test_Inj2(1:2:end)) + abs(Test_Inj2(2:2:end));
plot(3,abs(Jv(3)),'*','markersize',6,'color','red','LineWidth',0.6)
hold on
plot([1:2 4:179],abs(Jv([1:2 4:179])),'*','markersize',4,'color','black','LineWidth',0.6)
hold on
for ii = 1:179
    line([ii; ii],[1e-15; Jv(ii)],'Color','black')
end
set(gcf,'Units','inches','Position',[0 0 9 2.5])
set(gca,'FontName','Times','FontSize',14)
ylim([0 0.15])
xlabel({'${\rm Generator \; Index}$'},'Interpreter','latex','FontSize',14)
ylabel({'${\rm Current \; Injection\;}\overline{\bf J}_i$'},'Interpreter','latex','FontSize',14)
text(6,0.1328,'$\leftarrow {\rm FO \; Current \; Injection}$','Interpreter','latex','FontSize',14)
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"
