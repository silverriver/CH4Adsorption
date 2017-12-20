function [r, nex_calculated, nab_calculated, nbulk_calculated, sigma, a, output] = cal_r (nex_experiment, rou_free)
%nex_experiment (Unit:mmol/g) is the adsorption amount obtained using the experiment
%rou_free       (Unit:g/cm3)  is the density of gas at the free state

%check if the input have the same length
if length(nex_experiment) ~= length(rou_free)
    disp nex and rou_free do not have the same length
    return 
end

%calculate the value of r using the optimization algorithms
r = zeros(1,length(nex_experiment));
for i=1:length(nex_experiment)
    r(i) = fminsearch (@(x) abs(cal_adsorption(x, rou_free(i))-nex_experiment(i)), 10^(-6), optimset('Tolfun', 1e-8));
end

%calculate other values using r
nex_calculated   = zeros(1,length(nex_experiment));
nab_calculated   = zeros(1,length(nex_experiment));
nbulk_calculated = zeros(1,length(nex_experiment));
sigma            = zeros(1,length(nex_experiment));
a                = zeros(1,length(nex_experiment));
for i=1:length(nex_experiment)
    [nex_calculated(i), nab_calculated(i), nbulk_calculated(i), sigma(i), a(i)] = cal_adsorption(r(i),rou_free(i));
end
output = [nab_calculated; nbulk_calculated; nex_calculated; r; a; sigma];

end 