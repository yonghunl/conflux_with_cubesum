close all; clear all; clc;
%% Liveness Attack Scenario for reation time 1
% honest nodes - following the heaviest chain
% adv. nodes - keeping the balance b/t the 1st & 2nd heaviest chain
% Compute and compare the reation time 1 for conflux and cubesum 
% according to varying the adv. ratio

chain_A = [];
chain_B = [];

time_stamp = 0;

%% Simulation Paremeter setting

% Cubesum parameters
% Delta_threshold = 30;
% B_threshold = 10;
Delta_threshold = 1800; % corresponding to alpha
B_threshold = 400;


% Conflux parameters
h= 600;
h_0 = 360;
alpha = 1800;
beta = 160;
gamma = 10000;

N_total_blc = 15000;
ADV = -1;
HON = 1;

% Block generation rate
blc_gen_rate = 4; % (honest + adv) gen rate
adv_ratio_arr = [0.1, 0.15, 0.2, 0.25, .3, .35, 0.4, 0.5];
% adv_ratio_arr = [.35, .4, .45, .5];

avg_reaction_1 = zeros(2, length(adv_ratio_arr));

for n_adv_ratio = 1: length(adv_ratio_arr)
%     n_adv_ratio = 1;
%     adv_ratio = 0.5; % begin with the simplest case
    adv_ratio = adv_ratio_arr(n_adv_ratio); % begin with the simplest case
    
    adv_gen_rate = blc_gen_rate * adv_ratio;
    hon_gen_rate = blc_gen_rate * (1-adv_ratio);
        
    N_rep_adv = 50;
    reaction_time_1 = zeros(2, N_rep_adv); % 1st row for cubesum; 2nd row for conflux
    
    for n_rep_each_adv_ratio = 1: N_rep_adv

        % exponential distibution is memoryless 
        adv_blc_gen_time = exprnd(1/adv_gen_rate, 1, N_total_blc);
        hon_blc_gen_time = exprnd(1/hon_gen_rate, 1, N_total_blc);

        Delta_cube = 0;
        B_cube = 0;
        current_flag_cube = 0; % 1: under liveness attack, 0: normal situation
        current_flag_conflux = 0; % 1: under liveness attack, 0: normal situation

        Cubesum_array = [2, N_total_blc]; %[Delta; B]
        Conflux_array = [3, N_total_blc]; %[f(a); g(a); t(a)_]

        % Reset the parameters
        Time_M1_to_M2_cube = 0;
        Time_M1_to_M2_conflux = 0;
        chain_A = struct([]);
        chain_B = struct([]);
        time_stamp = 0;
        
        for n_blc = 1: N_total_blc

            % Generate a new block
            time_stamp = time_stamp + min(adv_blc_gen_time(n_blc), hon_blc_gen_time(n_blc));
            new_blc.time_stamp = time_stamp; 

            new_blc.weight = 1; %(1:honest , -1: adv)
            if adv_blc_gen_time(n_blc) < hon_blc_gen_time(n_blc)
                new_blc.type = ADV; %(1:honest , -1: adv)
            else 
                new_blc.type = HON; %(1:honest , -1: adv)    
            end
            
            % assign the default value to flags
            new_blc.flag_cube = 0;
            new_blc.flag_conflux = 0;
            
            % calculate weight_A weight_B when there is no adaptive weights
            weight_A_prev = length(chain_A);
            weight_B_prev = length(chain_B);
            
            % Decide which chain the honest/adv. choose
            [chain_A, chain_B, weight_A_new, weight_B_new] = choose_chain(new_blc, chain_A, chain_B, weight_A_prev, weight_B_prev);

            % Cube_sum - time to notice a liveness attack
            [new_blc.flag_cube, Delta_cube, B_cube] = cube_sum_check(chain_A, chain_B, B_cube, Delta_threshold, B_threshold, weight_A_new, weight_B_new);
            Cubesum_array(1,n_blc) = Delta_cube;
            Cubesum_array(2,n_blc) = B_cube;
            if (n_blc >1) && (Time_M1_to_M2_cube == 0)
                if (new_blc.flag_cube == 1) && (prev_blc.flag_cube == 0)
                    Time_M1_to_M2_cube = new_blc.time_stamp;
                end
            end

            % Conflux - time to notice a liveness attack
            new_blc.flag_conflux = conflux_check(chain_A, chain_B, alpha, beta, gamma, new_blc.time_stamp, weight_A_new, weight_B_new);
            if n_blc >1 && (Time_M1_to_M2_conflux == 0)
                if (new_blc.flag_conflux == 1) && (prev_blc.flag_conflux == 0)
                     Time_M1_to_M2_conflux = new_blc.time_stamp;
                end
            end
            
            if (Time_M1_to_M2_cube > 0) && (Time_M1_to_M2_conflux > 0 )
                reaction_time_1(1,n_rep_each_adv_ratio) = Time_M1_to_M2_cube;
                reaction_time_1(2,n_rep_each_adv_ratio) = Time_M1_to_M2_conflux;
                break; 
            end
            
            prev_blc = new_blc;
            
        end % for of n_blc
        
    end % for of n_rep

    avg_reaction_time_1_cube = mean(reaction_time_1(1,:));
    avg_reaction_time_1_conflux = mean(reaction_time_1(2,:));
    
    avg_reaction_1(1,n_adv_ratio) = avg_reaction_time_1_cube;
    avg_reaction_1(2,n_adv_ratio) = avg_reaction_time_1_conflux;
end













