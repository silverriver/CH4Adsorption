function [nex, nab, nbulk, sigma, a] = cal_adsorption (r, rou_free)
%calculate the adsorption amount if r and rou are given
%r         (Unit:cm)    is the thickness of the adsorption layer
%rou_free  (Unit:g/cm3) is the density of gas in free state

%the returning values:
%nex       (Unit:mmol/g)
%nab       (Unit:mmol/g)
%nbulk     (Unit:mmol/g)
%sigma     (No Unit) parameter for the distribution function
%a         (No Unit) parameter for the distribution function

mu  = 0.58*10^(-7);             %(Unit: cm)    mu is a constant value. (mu is used as a constant value in the normal distribution)
rou_min = rou_free;            %(Unit: g/cm3) the density of gas in free state
rou_max = 0.421;               %(Unit: g/cm3) the maximum density of adsorbed gas
s = 6.24*10^(4);                %(Unit: cm2/g) surface area per unit mass

sigma = r/sqrt(2*log(rou_max/rou_min));     %(sigma is used as a constant value in the normal distribution)
a = rou_max*sigma*sqrt(2*pi);               %a is the constant to adjust the value of the normal distribution function

nbulk = s*r*rou_min;
nab   = s*a*(normcdf(r+mu,mu,sigma) - 0.5);
nex   = nab-nbulk;                          

%convert the unit from (g/g) to (mmol/g)
nex = nex*1000/16;
nab = nab*1000/16;
nbulk = nbulk*1000/16;
end