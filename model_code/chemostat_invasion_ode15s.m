function [current_ET,par_store, opt_growth,total_rho,species_enzymes,simu,num_auto_peaks] ...
    = chemostat_invasion_ode15s()

%This script runs the simulation of the evolutionary chemostat model

global par

%Set-up enzyme storage
species_enzymes = [];

%Set-up system with fresh or existing environment
if ~par.existing_environment
    disp('Starting invasions in pristine environment')
    simu.rho = [];
    simu.c = zeros(par.m+1,par.n_rxn+1);
    simu.c(end,1) = par.supply/par.D;
    if isfield(par,'supply2')
        simu.c(end,2) = simu.c(end,2) + par.supply2;
    end
    
else
    disp('Starting invasions in existing environment')
    simu.rho = par.existing_rho;
    simu.c = par.existing_c;
end

if ~par.existing_environment
    %Add in first strain
    [optimal_ET,opt_init_growth,c_ss] = compute_optimal_strain(simu);
    disp(['Novel strain class invading with growth rate ',num2str(-opt_init_growth),'.']);
    simu.c(end+1,:) = simu.c(end,:);
    simu.c(end-1,:) = c_ss;
    par.E = [par.E; optimal_ET(1:par.n_rxn)];
    par.T = [par.T; optimal_ET(par.n_rxn+1:end)];
    species_enzymes = [species_enzymes; optimal_ET];
    simu.rho(end+1,1) = par.invasion_biomass;
    par.m = par.m + 1;
end

o = 0;
round = 0;
num_auto_peaks = [];
species_list = (1:size(par.E,1))';
species_counter = size(par.E,1);

%Make optimal growth cell
opt_growth = [];

while o == 0
    round = round + 1;
    par_store{round} = par;
    disp(['Round ', num2str(round),': now running simulation with ',num2str(par.m), ' species.']);
    [~,rho_final,c_store{round},rho_store{round},t_store{round},simu] = ...
        simulate_crossfeeding_ode15s(simu);
    species_store{round} = species_list;
    
    %Remove failed strains
    extinct_index = rho_final < par.extol;
    simu.rho = simu.rho(~extinct_index);
    par.E = par.E(~extinct_index,:);
    par.T = par.T(~extinct_index,:);
    simu.c = [simu.c(~extinct_index,:);simu.c(end,:)];
    par.m = par.m - sum(extinct_index);
    species_list = species_list(~extinct_index);
    disp([num2str(sum(extinct_index)),' species have gone extinct.'])
    
    %Check if simulation exceeded time and thus failed to reach SS
    if par.exceeded_time
        
        o = 1;
        current_ET = [par.E,par.T];
        disp('Failed to reach steady state. Simulation terminating.')

        %Compute autocorrelation of final half of timecourse or last 20,000
        %points, which is shorter
        max_timepoints = 2e5;
        max_time = t_store{round}(end);
        [~,mid_ind] = min(abs(t_store{round} - max_time/2));
        if length(t_store{round}) - mid_ind > max_timepoints
            mid_ind = length(t_store{round}) - max_timepoints;
        end
        f = rho_store{round}(1,mid_ind:end);
        t = t_store{round}(mid_ind:end);
        C = spline_autocorr(t,f);
        num_auto_peaks =  length(findpeaks(C));

    %Continue simulation if integrator reached steady state
    else
        
        %Add in new optimal strain if growth exceeds death rate
        if sum(extinct_index) < length(rho_final)
            [optimal_ET,opt_growth(round),c_ss] = compute_optimal_strain(simu);
            
            %Determine whether new strain is same class as old strain
            new_strain_type = optimal_ET > 1e-15;
            old_strain_type = [par.E,par.T] > 1e-15;
            
            [~,~,identical_type] = intersect(new_strain_type,old_strain_type,'rows');
            
            %Check if there is substantial growth rate improvement
            if round > 1
                invasion_delta = abs((opt_growth(round) - opt_growth(round - 1))/opt_growth(round-1));
                if invasion_delta < 1e-4
                    proceed_invasion = 0;
                end
            else
                proceed_invasion = 1;
            end
            
            if -opt_growth(round) > par.D*(1+par.growth_tol)
                
                %Add new strain class with set invasion biomass 
                if isempty(identical_type) && proceed_invasion
                    simu.c(end+1,:) = simu.c(end,:);
                    simu.c(end-1,:) = c_ss;
                    par.E = [par.E; optimal_ET(1:par.n_rxn)];
                    par.T = [par.T; optimal_ET(par.n_rxn+1:end)];
                    simu.rho(end+1,1) = par.invasion_biomass;
                    par.m = par.m + 1;
                    species_list = [species_list;species_counter + 1];
                    species_counter = species_counter + 1;
                    disp(['Novel strain class invading with growth rate ',num2str(-opt_growth(round)),'.']);
                
                %Add existing strain class by replacing previous strain of
                %that type, leaving biomass the same
                elseif ~isempty(identical_type) && proceed_invasion
                    simu.c(identical_type,:) = simu.c(end,:);
                    simu.c(identical_type,:) = c_ss;
                    par.E(identical_type,:) = optimal_ET(1:par.n_rxn);
                    par.T(identical_type,:) = optimal_ET(par.n_rxn+1:end);
                    %simu.rho(identical_type,1) = par.invasion_biomass;
                    disp(['Existing strain class invading with growth rate ',num2str(-opt_growth(round)),'.']);
                    
                    %End if no growth rate improvement
                else
                    o = 1;
                    current_ET = [par.E,par.T];
                    disp(['No growth rate improvement (',num2str(-opt_growth(round)),'). Simulation terminating.'])
                end
                
                species_enzymes = [species_enzymes; optimal_ET];
                
                %End if none can invade
            else
                o = 1;
                current_ET = [par.E,par.T];
                disp('No new species can invade. Simulation terminating')
            end
            
            %Exit if final round
            if round == par.n_rounds
                o = 1;
                current_ET = [par.E,par.T];
                disp('Final round completed. Simulation terminating.')
            end
            
        else
            o = 1;
            current_ET = [par.E,par.T];
            opt_growth = 0;
            disp('All species have gone extinct. Simulation terminating')
        end
    end
end



%Compute aggregate timecourses
total_time = sum(cellfun(@(x) size(x,2),rho_store));
total_c = [];
total_rho = zeros(species_counter,total_time);
total_t = zeros(1,total_time);
start_ind = 1;
for i = 1:round
    final_ind = size(rho_store{i},2);
    total_c = [total_c , squeeze(c_store{i}(end,:,:))];
    if i == 1
        total_t(start_ind:(final_ind+start_ind-1)) = t_store{i};
    else
        total_t(start_ind:(final_ind+start_ind-1)) = t_store{i}+total_t(start_ind-1);
    end
    species_i = species_store{i};
    for j = 1:length(species_i)
        total_rho(species_i(j),start_ind:(final_ind+start_ind-1)) = rho_store{i}(j,:);
    end
    start_ind = final_ind+start_ind;
end

if par.plt_invasion
    
    figure
    plot(total_t,total_rho','LineWidth',1.5)
    biomass_leg_entries = cellstr(cellfun(@(x) num2str(x),num2cell((1:species_counter)'),'UniformOutput',false));
    biomass_leg_entries = strcat('$\rho_{',biomass_leg_entries);
    biomass_leg_entries = strcat(biomass_leg_entries,'}$');
    legend(biomass_leg_entries,'Interpreter','latex')
    title('Cell counts','Interpreter','latex')
    ylabel('Cell counts','Interpreter','latex')
    xlabel('Time','Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    %set(gca,'YScale','log')
    if par.save_invasion
        print(gcf, '-dpng',['plots/invasion_biomass_conc.png'],'-r400');
    end
    
    
    figure
    biomass_leg_entries2 = cellstr(cellfun(@(x) num2str(x),num2cell((1:(round+1))'),'UniformOutput',false));
    biomass_leg_entries2 = strcat('$\rho_{',biomass_leg_entries2);
    biomass_leg_entries2 = strcat(biomass_leg_entries2,'}$');
    imagesc(species_enzymes)
    xticks(1:size(species_enzymes,2))
    xticklabels({'$E_1$','$E_2$','$T_1$','$T_2$','$T_3$'})
    yticks(1:size(species_enzymes,1))
    yticklabels(biomass_leg_entries2)
    set(gca,'TickLabelInterpreter','latex')
    title('Species strategies','Interpreter','latex')
    set(gca,'FontSize',14)
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    if par.save_invasion
        print(gcf, '-dpng',['plots/species_strategies.png'],'-r400');
    end
end

end

