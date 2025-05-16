function [c_final,rho_final,c_store,rho_store,t_store,simu] = simulate_crossfeeding_ode15s(simu)
global par

%This script runs a purely ecological chemostat simulation, it is called by
%"chemostat_invasion_ode15s.m"

%Apply defaults if not specified
if ~isfield(par,'no_saturation')
    par.no_saturation = false;
end
if ~isfield(par,'no_thermo')
    par.no_thermo = false;
end

if ~isfield(par,'max_internal_time')
    par.max_internal_time = 1e13;
end

y0 = [simu.rho;simu.c(:)];
y0(y0<0) = 0;
opt_cell{1} = odeset('Events', @SS_criterion,'NonNegative',1:length(y0),'RelTol',1e-9,'AbsTol',1e-6);
opt_cell{2} = odeset('Events', @SS_criterion,'NonNegative',1:length(y0));
tspan = [0,par.max_internal_time];
odefun = @(t,y) crossfeeding_RHS_ode15s(y);

% par.start_time = datetime('now');
% [t_store,y] = ode15s(odefun,tspan,y0,opt_cell{1});
% %[t_store,y] = ode45(odefun,tspan,y0,Opt);

%Rerun with different integrator settings if steady-state convergence not reached
loop_var = true;
tn = 1;
while loop_var

    %Get starting time for timeout and run
     par.start_time = datetime('now');
     par.exceeded_time = false;
     [t_store,y] = ode15s(odefun,tspan,y0,opt_cell{tn});

     %End if internal time not exceeded or final integrator setting reached
     if ~par.exceeded_time || tn == length(opt_cell)
         loop_var = false;
     %End if internal time exceeded and final integrator
     elseif par.exceeded_time && tn == length(opt_cell)
         par.exceeded_runtime = 1;
         loop_var = false;
     end

     tn = tn + 1;
end


%Arrange storage variables
rho_store = y(:,1:par.m)';
rho_final = rho_store(:,end);
simu.rho = rho_final;

c_store = zeros(par.m+1,par.n_rxn+1,length(t_store));
c_data = y(:,par.m+1:end);
for i = 1:length(t_store)
    c_store(:,:,i) = reshape(c_data(i,:),par.m+1,par.n_rxn+1);
end
c_final = c_store(:,:,end);
simu.c = c_final;

if par.plt
    
    leg_entries = cellstr(cellfun(@(x) num2str(x),num2cell((1:size(c_final,2))')));
    leg_entries = strcat('$c_',leg_entries);
    leg_entries = strcat(leg_entries,'$');
    for i = 1:(par.m+1)
        figure
        data = squeeze(c_store(i,:,:))';
        plot(t_store,data,'LineWidth',1.5)
        
        if i == (par.m+1)
            title('Extracellular concentrations','Interpreter','latex')
        else
            title(['Intracellular concentrations of species ',num2str(i)],'Interpreter','latex')
        end
        legend(leg_entries,'Interpreter','latex')
        ylabel('Concentration','Interpreter','latex')
        xlabel('Time','Interpreter','latex')
        set(gca,'FontSize',14)
        set(gca,'TickLabelInterpreter','latex')
        
        if par.save
            if i == (par.m+1)
                print(gcf, '-dpng',['plots/extracellular_conc.png'],'-r400');
            else
                print(gcf, '-dpng',['plots/species_',num2str(i),'_conc.png'],'-r400');
            end
        end
    end
    
    figure
    plot(t_store,rho_store','LineWidth',1.5)
    biomass_leg_entries = cellstr(cellfun(@(x) num2str(x),num2cell((1:size(c_final,1)-1)')));
    biomass_leg_entries = strcat('$\rho_',biomass_leg_entries);
    biomass_leg_entries = strcat(biomass_leg_entries,'$');
    legend(biomass_leg_entries,'Interpreter','latex')
    title('Cell counts','Interpreter','latex')
    ylabel('Concentration','Interpreter','latex')
    xlabel('Time','Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'YScale','log')
    if par.save
        print(gcf, '-dpng',['plots/biomass_conc.png'],'-r400');
    end
end



end

