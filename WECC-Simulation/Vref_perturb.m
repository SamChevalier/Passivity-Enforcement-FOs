function Vref_perturb(t)
global Source_V Amplitude_V Frequency_V Exc vref0

global intstep VP_nos VQ_nos tcorr P0 Q0 PQ_std Fl ij
   
% 1. Cause the Reference to Oscillate on Bus 4 (Index 3)
    Exc.vref0(Source_V) = vref0*(1 + (Amplitude_V*(1+0.025*randn))*sin(2*pi*Frequency_V*t + 0.025*randn));

% 2. Add Load Perturbations
    gamma    = 1/tcorr;                  % Inverse noise correlation time
    nos_std  = PQ_std*sqrt(2*gamma);     % Should be 0.01*sqrt(2*gamma); 
    rnd_vec  = randn(length(P0),1);
    
    % Load noise
    VP_nos(:,ij+1) = VP_nos(:,ij)*(1 - gamma * intstep) + nos_std * sqrt(intstep) * rnd_vec;
    VQ_nos(:,ij+1) = VQ_nos(:,ij)*(1 - gamma * intstep) + nos_std * sqrt(intstep) * rnd_vec;
    ij             = ij + 1;
    
    % Add the noise
    Fl.con(:,2) = P0.*(1+VP_nos(:,ij));
    Fl.con(:,5) = Q0.*(1+VQ_nos(:,ij));
end