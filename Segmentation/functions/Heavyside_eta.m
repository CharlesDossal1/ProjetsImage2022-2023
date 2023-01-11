%Approximation of the dirac function
function [M] = Heavyside_eta(phi,eta)
M = (1.+2.*atan(phi./eta)/pi)./2.;