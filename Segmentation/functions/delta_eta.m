%Approximation of the dirac function
function [M] = delta_eta(phi,eta)
M = 1./((1 + (phi/eta).^2)*eta*pi);