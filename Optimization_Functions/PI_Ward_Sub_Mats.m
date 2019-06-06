function [m1,m2,m3,m4,Ybn,bus_vec] = PI_Ward_Sub_Mats(Yb,inj_bus,edge_bus)
% PI_WARD_SUB_MATS   Construct "m" matrices
%
% Inputs:   
% 1) Yb              Full (2n)x(2n) Y-bus matrix, where the shunt elements
%                    at bus ind_inj and bus(es) ind_edge have been removed
% 2) ind_inj         Bus to perform the reduction at
% 3) ind_edge        Buses whose injections are given
%
% Outputs:
% 1) Yeq             Dynamic Ward Equivalent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(Yb)/2;
m = length(edge_bus);
bus_vec = 1:n;
bus_vec([inj_bus edge_bus]) = [];

% Create vector ordered correctly
bus_vec = [inj_bus edge_bus bus_vec];

% Convert to augmented form
bus_vec_aug          = zeros(1,2*length(bus_vec));
bus_vec_aug(1:2:end) = 2*bus_vec-1;
bus_vec_aug(2:2:end) = 2*bus_vec;

% Reorder the Yb matrix
Ybn = Yb;
Ybn = Ybn(bus_vec_aug,:);
Ybn = Ybn(:,bus_vec_aug);

% Define sub-Y matrcies
Y1 = Ybn(1:2         ,1:2);
Y2 = Ybn(1:2         ,3:(2*m+2));
Y3 = Ybn(1:2         ,(2*m+3):end);
Y4 = Ybn(3:(2*m+2)   ,1:2);
Y5 = Ybn(3:(2*m+2)   ,3:(2*m+2));
Y6 = Ybn(3:(2*m+2)   ,(2*m+3):end);
Y7 = Ybn((2*m+3):end ,1:2);
Y8 = Ybn((2*m+3):end ,3:(2*m+2));
Y9 = Ybn((2*m+3):end ,(2*m+3):end);

% Define m matrices
m1 = Y1 - Y3*(Y9\Y7);
m2 = Y2 - Y3*(Y9\Y8);
m3 = Y4 - Y6*(Y9\Y7);
m4 = Y5 - Y6*(Y9\Y8);
end

