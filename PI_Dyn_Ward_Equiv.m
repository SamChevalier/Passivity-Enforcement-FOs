function [Yeq] = PI_Dyn_Ward_Equiv(Yb,bus)
% PI_DYN_WARD_EQUIV  Build Dynamic Ward Equivalent
%
% Inputs:   
% 1) Yb              Full (2n)x(2n) Y-bus matrix (Shunt at source excluded)
% 2) bus             Bus to perform the reduction at
%
% Outputs:
% 1) Yeq             Dynamic Ward Equivalent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(Yb)/2;

if bus == 1
    Y1  = Yb(1:2,1:2);
    Y2  = Yb(1:2,3:end);
    Y3  = Yb(3:end,1:2);
    Y4  = Yb(3:end,3:end);
    Yeq = Y1-Y2*(Y4\Y3);
    
elseif bus == n
    Y1  = Yb(1:(end-2),1:(end-2));
    Y2  = Yb(1:(end-2),(end-1):end);
    Y3  = Yb((end-1):end,1:(end-2));
    Y4  = Yb((end-1):end,(end-1):end);
    Yeq = Y4-Y3*(Y1\Y2);
    
else
    Y1 = Yb(1:(2*bus-2)       ,1:(2*bus-2));
    Y2 = Yb(1:(2*bus-2)       ,(2*bus-1):(2*bus));
    Y3 = Yb(1:(2*bus-2)       ,(2*bus+1):end);
    Y4 = Yb((2*bus-1):(2*bus) ,1:(2*bus-2));
    Y5 = Yb((2*bus-1):(2*bus) ,(2*bus-1):(2*bus));
    Y6 = Yb((2*bus-1):(2*bus) ,(2*bus+1):end);
    Y7 = Yb((2*bus+1):end     ,1:(2*bus-2));
    Y8 = Yb((2*bus+1):end     ,(2*bus-1):(2*bus));
    Y9 = Yb((2*bus+1):end     ,(2*bus+1):end);
    
    % Build new matrix
    Ybn = [Y5 Y4 Y6;
           Y2 Y1 Y3;
           Y8 Y7 Y9];
    
    % Perform standard Kron reduction
    Y1  = Ybn(1:2,1:2);
    Y2  = Ybn(1:2,3:end);
    Y3  = Ybn(3:end,1:2);
    Y4  = Ybn(3:end,3:end);
    Yeq = Y1-Y2*(Y4\Y3);
end

