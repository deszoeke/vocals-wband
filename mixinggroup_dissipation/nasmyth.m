function [phi11,phi22] = nasmyth(k,epsilon,nu)

% NASMYTH Nasmyth universal velocity spectrum.
%
% Inputs
% k --> wavenumber in cpm
% epsilon --> dissipation
% nu --> viscosity of seawater (default 1e-6 m^2/s)

% Outputs
% phi11 --> Spectrum for U1 (vel. in dir. of k)
% phi22 --> Spectrum for U2 (vel. perpindicular to k)
%
% The form of the spectrum is computed from a spline interpolation of 
% Nasmyth points listed by Oakey (1982).
% 
% References: Oakey, N. S., 1982: J. Phys. Ocean., 12, 256-271.
         

if nargin<3
    nu=1e-6;
end

kks=[2.83e-4 5.03e-4 8.95e-4 1.59e-3 2.83e-3 5.03e-3 8.95e-3 ...
    1.59e-2 2.83e-2 5.03e-2 7.97e-2 1.26e-1 1.59e-1 2.00e-1 2.52e-1];
F1=[1.254e5 4.799e4 1.842e4 7.050e3 2.699e3 1.036e3 3.964e2 ...
    1.490e2 3.574e1 5.600e0 7.214e-1 6.580e-2 1.812e-2 4.552e-3 1.197e-3];
F2=[1.671e5 6.397e4 2.455e4 9.404e3 3.598e3 1.380e3 5.320e2 ...
    2.302e2 6.880e1 1.373e1 2.101e0 2.127e-1 6.161e-2 1.570e-2 4.011e-3];

ks=(epsilon./nu.^3).^(1/4);

F1spline=10.^spline(log10(kks),log10(F1),log10(k./ks));
F2spline=10.^spline(log10(kks),log10(F2),log10(k./ks));
phi11=(epsilon*nu^5).^0.25.*F1spline;
phi22=(epsilon*nu^5).^0.25.*F2spline;
