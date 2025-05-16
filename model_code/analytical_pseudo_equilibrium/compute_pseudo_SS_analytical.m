function [pseudo_growth,growth,c_ss] = compute_pseudo_SS_analytical(ET,simu)

%Compute the pseudo-equilbrium of the model with linear kinetics and
%thermodynamic inhibition. This version does not consider the initial
%effect of growth dilution

global par

%Determine which solution to use based on enzyme allocation
if ET(2) == 0 && ET(5) == 0
    c_ss = ...
        P1_analytical_pseudo_equi(ET(1),ET(2),ET(3),ET(4),ET(5),...
        simu.c(end,1),simu.c(end,2),simu.c(end,3),par.nu,...
        par.beta,par.Keq(1),par.Keq(2));
elseif ET(4) == 0
    c_ss = ...
        P2_analytical_pseudo_equi(ET(1),ET(2),ET(3),ET(4),ET(5),...
        simu.c(end,1),simu.c(end,2),simu.c(end,3),par.nu,...
        par.beta,par.Keq(1),par.Keq(2));
elseif ET(1) == 0 && ET(3) == 0
    c_ss = ...
        C_analytical_pseudo_equi(ET(1),ET(2),ET(3),ET(4),ET(5),...
        simu.c(end,1),simu.c(end,2),simu.c(end,3),par.nu,...
        par.beta,par.Keq(1),par.Keq(2));
else
    c_ss = ...
        vers_analytical_pseudo_equi(ET(1),ET(2),ET(3),ET(4),ET(5),...
        simu.c(end,1),simu.c(end,2),simu.c(end,3),par.nu,...
        par.beta,par.Keq(1),par.Keq(2));
end

%Compute the growth rates
Jg1 = max([ET(1)*(par.nu)*(c_ss(1) - c_ss(2)/par.Keq(1)),0]);
Jg2 = max([ET(2)*(par.nu)*(c_ss(2) - c_ss(3)/par.Keq(2)),0]);
total_Jg = max([par.alpha(1)*Jg1 + par.alpha(2)*Jg2 - par.min_growth,0]);
pseudo_growth = (1 - par.osm_coeff*(c_ss(1) + c_ss(2) +c_ss(3)))*total_Jg;
pseudo_growth(isnan(pseudo_growth)) = -100;
growth = max([pseudo_growth,0]);

end

