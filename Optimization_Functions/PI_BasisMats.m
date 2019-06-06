function [BasisCoeffs] = PI_BasisMats(Y0)
% PI_BASISMATS:   Take numerical matrix Y0 and deconstruct it into its 8
%                 component basis matrices.
%
% Inputs:         
% 1) Y0           Numerical admittance matrix
%
% Outputs:
% 1) BasisCoeffs  A structure of coefficients needed to reconstruct Y0, 
%                 such that Y0 = sum ((ai + j*bi)*Ti):
%                     1) BasisCoeffs.aT1
%                     2) BasisCoeffs.aT2
%                     3) BasisCoeffs.aT3
%                     4) BasisCoeffs.aT4
%                     5) BasisCoeffs.bT1
%                     6) BasisCoeffs.bT2
%                     7) BasisCoeffs.bT3
%                     8) BasisCoeffs.bT4
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct 8 linear equations and with 8 unknowns: Ax = b, where we have
% that x = [aT1; aT2; aT3; aT4; bT1; bT2; bT3; bT4];
b = [real(Y0(1,1));
     imag(Y0(1,1));
     real(Y0(1,2));
     imag(Y0(1,2));
     real(Y0(2,1));
     imag(Y0(2,1));
     real(Y0(2,2));
     imag(Y0(2,2))];

% Static "A" matrix
A = [1 1 0  0 0 0 0 0;
     0 0 0  0 1 1 0 0;
     0 0 1 -1 0 0 0 0;
     0 0 0 0 0 0 1 -1;
     0 0 1 1 0 0 0  0;
     0 0 0 0 0 0 1  1;
     1 -1 0 0 0 0 0 0;
     0 0 0 0 1 -1 0 0];
 
% Solve the system
x = A\b;

% Parse
BasisCoeffs.aT1 = x(1);
BasisCoeffs.aT2 = x(2);
BasisCoeffs.aT3 = x(3);
BasisCoeffs.aT4 = x(4);
BasisCoeffs.bT1 = x(5);
BasisCoeffs.bT2 = x(6);
BasisCoeffs.bT3 = x(7);
BasisCoeffs.bT4 = x(8);

end

