function S = source(S0, u)
%This function defines the relation between the fuel concentration f and 
%the source term S(f)
%INPUTS:
%  S0 = reaction rate amplitude
%  u = state vectors for all elements (rho, rho*u, rho*v, rho*E, rho*f)
%OUTPUT:
%  S = reaction rate

h = 10; A = 100; b = 0.5; fMax = 1;

f = u(5) / u(1);
if f > 0 && f < fMax
    Sf = S0 * f * (fMax - f) * exp(-A * (f / fMax - b)^2);
else
    Sf = 0;
end
S = [0, 0, 0, -h * Sf, Sf];

end