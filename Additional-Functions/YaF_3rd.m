function [YaF] = YaF_3rd(~)
% YAF_3RD_AVR     Build the FRF associated with a 3rd order generator
%
% Inputs:        NA
%
% Outputs:
% 1) YaF         FRF in functional form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j = sqrt(-1);

% Define Variables
syms Vr Vi Omega Pm Xd Xdp Xq Td0p H D Ef
syms d w eqp

% Voltage Magnitude
V = sqrt(Vr^2 + Vi^2);

% Voltages and Currents
ed = sin(d)*Vr - cos(d)*Vi;
eq = cos(d)*Vr + sin(d)*Vi;
id = (eqp - eq)/Xdp;
iq = ed/Xq;

% Outputs
Ir =  sin(d)*id + cos(d)*iq;
Ii = -cos(d)*id + sin(d)*iq;

% Electrical Powers
Pe = ed*id+eq*iq;

% Vectors
X  = [d; w; eqp];
U  = [Vr; Vi];
G  = [Ir; Ii];

% System DAE Model
w0 = 2*pi*60;
F = [w;
    (Pm - Pe - D*w)/(2*H/w0);
    (Ef-(Xd-Xdp)*id-eqp)/Td0p];

% Build Jacobians
JFX = jacobian(F,X);
JFU = jacobian(F,U);
JGX = jacobian(G,X);
JGU = jacobian(G,U);

% Output Matrix
Ya  = ((JGX/(j*Omega*eye(length(X)) - JFX))*JFU + JGU);
YaF = matlabFunction(Ya);

end

