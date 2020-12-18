close all; clear all; clc;
%% Liveness Attack Scenario for reation time
% honest nodes - following the heaviest chain
% adv. nodes - keeping the balance b/t the 1st & 2nd heaviest chain
% Compute and compare the reation time 1 for conflux and cubesum 
% according to varying the adv. ratio

chain_A = [];
chain_B = [];

time_stamp = 0;

%% Simulation Paremeter setting

% Cubesum parameters
Delta_threshold = 1800; % corresponding to alpha
B_threshold = 160;

% Conflux parameters
h= 600;
h_0 = 360;
alpha = 1800;
beta = 160;
gamma = 10000;

N_total_blc = 35000;
ADV = -1;
HON = 1;

CHAIN_A = 1;
CHAIN_B = 2;


% Block generation rate
blc_gen_rate = 4; % (honest + adv) gen rate
adv_ratio_arr = [0.3];

avg_reaction_1 = zeros(2, length(adv_ratio_arr));
avg_reaction_2 = zeros(2, length(adv_ratio_arr));

for n_adv_ratio = 1: length(adv_ratio_arr)

    % set the block generation ratio for hon./adv. 
    % according to the adv ratio
    adv_ratio  = adv_ratio_arr(n_adv_ratio);
    adv_gen_rate = blc_gen_rate * adv_ratio;
    hon_gen_rate = blc_gen_rate * (1-adv_ratio);

    reaction_time_1 = zeros(2, N_rep_adv); % 1st row for cubesum; 2nd row for conflux
    reaction_time_2 = zeros(2, N_rep_adv); % 1st row for cubesum; 2nd row for conflux
    
    N_rep_adv = 100;
    for n_rep_each_adv_ratio = 1: N_rep_adv
        
        % exponential distibution is memoryless 
        adv_blc_gen_time = exprnd(1/adv_gen_rate, 1, N_total_blc);
        hon_blc_gen_time = exprnd(1/hon_gen_rate, 1, N_total_blc);

        Delta_cube = 0;
        B_cube = 0;
        current_flag_cube = 0; % 1: under liveness attack, 0: normal situation
        current_flag_conflux = 0; % 1: under liveness attack, 0: normal situation

        % cubesum reaction time 1 and 2
        Cubesum_array = [2, N_total_blc]; %[Delta; B]
        Time_M1_to_M2_cube = 0;
        Time_M2_to_M1_cube = 0;
        
        
        % Reset the parameters
        chain_A = struct([]);
        chain_B = struct([]);
        weight_A_prev = 0;
        weight_B_prev = 0;
        time_stamp = 0;
                
        prev_blc.flag_cube = 0;
        
        rand_for_h = rand(1,N_total_blc);
        cnt_adaptive_weight = 0;
        
        for n_blc = 1: N_total_blc % cubesum reaction time

            % Generate a new block
            time_stamp = time_stamp + min(adv_blc_gen_time(n_blc), hon_blc_gen_time(n_blc));
            new_blc.time_stamp = time_stamp; 
            
            if adv_blc_gen_time(n_blc) < hon_blc_gen_time(n_blc)
                new_blc.type = ADV; %(1:honest , -1: adv)
            else 
                new_blc.type = HON; %(1:honest , -1: adv)    
            end
            
            if new_blc.type == ADV
                new_blc.weight = 1; 
            else
                if prev_blc.flag_cube == 0 % normal weight
                    new_blc.weight = 1; 
                else % adaptive weigh
                    cnt_adaptive_weight = cnt_adaptive_weight + 1;
                    if rand_for_h(cnt_adaptive_weight) > 1/h
                        new_blc.weight = 0;
                    else
                        new_blc.weight = h; %(1:honest , -1: adv)
                    end
                end
            end

            % assign the default value to the flag
            new_blc.flag_cube = 0;
            
            % Decide which chain the honest/adv. choose
            if prev_blc.flag_cube == 0 %
                if Time_M1_to_M2_cube == 0
                    [chain_A, chain_B, weight_A_new, weight_B_new] = choose_chain_LA(new_blc, chain_A, chain_B, weight_A_prev, weight_B_prev);
                else
                    [chain_A, chain_B, weight_A_new, weight_B_new] = choose_chain(new_blc, chain_A, chain_B, weight_A_prev, weight_B_prev);
                end
            else
                if  new_blc.type == ADV % adv. blc
                    if chain_for_adv == CHAIN_A
                        chain_A = [chain_A, new_blc];
                        weight_A_new = weight_A_prev + new_blc.weight;
                    else
                        chain_B = [chain_B, new_blc];
                        weight_B_new = weight_B_prev + new_blc.weight;
                    end
                else % honest. blc
                    if chain_for_hon == CHAIN_A
                        chain_A = [chain_A, new_blc];
                        weight_A_new = weight_A_prev + new_blc.weight;
                    else
                        chain_B = [chain_B, new_blc];
                        weight_B_new = weight_B_prev + new_blc.weight;
                    end
                end
            end
            
            % Cube_sum - time to notice a liveness attack
            [new_blc.flag_cube, Delta_cube, B_cube] = cube_sum_check(chain_A, chain_B, B_cube, Delta_threshold, B_threshold, weight_A_new, weight_B_new);
            Cubesum_array(1,n_blc) = Delta_cube;
            Cubesum_array(2,n_blc) = B_cube;
            
            if n_blc >1
                if (new_blc.flag_cube == 1) && (prev_blc.flag_cube == 0) ...
                        && (Time_M1_to_M2_cube == 0) && (Time_M2_to_M1_cube == 0)
                    Time_M1_to_M2_cube = new_blc.time_stamp;
                    if weight_A_new >= weight_B_new
                        chain_for_hon = CHAIN_A;
                        chain_for_adv = CHAIN_B;
                    else
                        chain_for_hon = CHAIN_B;
                        chain_for_adv = CHAIN_A;
                    end
                end
                if (new_blc.flag_cube == 0) && (prev_blc.flag_cube == 1)...
                        && (Time_M1_to_M2_cube > 0) && (Time_M2_to_M1_cube == 0)
                    Time_M2_to_M1_cube = new_blc.time_stamp;
                end
            end

            if (Time_M1_to_M2_cube > 0) && (Time_M2_to_M1_cube > 0)
                reaction_time_1(1,n_rep_each_adv_ratio) = Time_M1_to_M2_cube;
                reaction_time_2(1,n_rep_each_adv_ratio) = Time_M2_to_M1_cube - Time_M1_to_M2_cube;
                break;
            end
            
            prev_blc = new_blc;
            weight_A_prev = weight_A_new; 
            weight_B_prev = weight_B_new; 
        end % cubesum
    

        % conflux reaction time 1 and 2
        Conflux_array = [3, N_total_blc]; %[f(a); g(a); t(a)]
        Time_M1_to_M2_conflux = 0;
        Time_M2_to_M1_conflux = 0;
        
        % Reset the parameters
        chain_A = struct([]);
        chain_B = struct([]);
        weight_A_prev = 0;
        weight_B_prev = 0;
        time_stamp = 0;
        
        prev_blc.flag_conflux = 0;
        cnt_adaptive_weight = 0;
        for n_blc = 1: N_total_blc % conflux reaction time

            % Generate a new block
            time_stamp = time_stamp + min(adv_blc_gen_time(n_blc), hon_blc_gen_time(n_blc));
            new_blc.time_stamp = time_stamp; 
            % new_blc.weight = 1; %(1:honest , -1: adv)
            
            if adv_blc_gen_time(n_blc) < hon_blc_gen_time(n_blc)
                new_blc.type = ADV; %(1:honest , -1: adv)
            else 
                new_blc.type = HON; %(1:honest , -1: adv)    
            end
            
            if  new_blc.type == ADV % normal weight
                new_blc.weight = 1; 
            else
                if prev_blc.flag_conflux == 0 % normal weight
                    new_blc.weight = 1; 
                else % adaptive weigh
                    cnt_adaptive_weight = cnt_adaptive_weight + 1;
                    if rand_for_h(cnt_adaptive_weight) > 1/h
                        new_blc.weight = 0;
                    else
                        new_blc.weight = h; %(1:honest , -1: adv)
                    end
                end
            end
            
            % assign the default value to the flag
            new_blc.flag_conflux = 0;
            
            % Decide which chain the honest/adv. choose
            if prev_blc.flag_conflux == 0 % normal mode
                if Time_M1_to_M2_conflux == 0
                    [chain_A, chain_B, weight_A_new, weight_B_new] = choose_chain_LA(new_blc, chain_A, chain_B, weight_A_prev, weight_B_prev);
                else
                    [chain_A, chain_B, weight_A_new, weight_B_new] = choose_chain(new_blc, chain_A, chain_B, weight_A_prev, weight_B_prev);
                end
            else % conservative mode (mode 2)
                if  new_blc.type == ADV % adv. blc
                    if chain_for_adv == CHAIN_A
                        chain_A = [chain_A, new_blc];
                        weight_A_new = weight_A_prev + new_blc.weight;
                    else
                        chain_B = [chain_B, new_blc];
                        weight_B_new = weight_B_prev + new_blc.weight;
                    end
                else % honest. blc
                    if chain_for_hon == CHAIN_A
                        chain_A = [chain_A, new_blc];
                        weight_A_new = weight_A_prev + new_blc.weight;
                    else
                        chain_B = [chain_B, new_blc];
                        weight_B_new = weight_B_prev + new_blc.weight;
                    end
                end
                    
            end
            
            % Conflux - time to notice a liveness attack
            new_blc.flag_conflux = conflux_check(chain_A, chain_B, alpha, beta, gamma, new_blc.time_stamp, weight_A_new, weight_B_new);
            
            if n_blc >1
                if (new_blc.flag_conflux == 1) && (prev_blc.flag_conflux == 0) ...
                        && (Time_M1_to_M2_conflux == 0) && (Time_M2_to_M1_conflux == 0)
                    Time_M1_to_M2_conflux = new_blc.time_stamp;
                    if weight_A_new >= weight_B_new
                        chain_for_hon = CHAIN_A;
                        chain_for_adv = CHAIN_B;
                    else
                        chain_for_hon = CHAIN_B;
                        chain_for_adv = CHAIN_A;
                    end
                    
                end
                if (new_blc.flag_conflux == 0) && (prev_blc.flag_conflux == 1) ...
                        && (Time_M2_to_M1_conflux == 0) && (Time_M1_to_M2_conflux > 0)
                    Time_M2_to_M1_conflux = new_blc.time_stamp;
                end
            end

            if (Time_M2_to_M1_conflux > 0) && (Time_M1_to_M2_conflux > 0)
                reaction_time_1(2,n_rep_each_adv_ratio) = Time_M1_to_M2_conflux;
                reaction_time_2(2,n_rep_each_adv_ratio) = Time_M2_to_M1_conflux - Time_M1_to_M2_conflux;
                break;
            end
                                   
            prev_blc = new_blc;
            weight_A_prev = weight_A_new; 
            weight_B_prev = weight_B_new; 
            
            if n_blc == N_total_blc - 1
                reaction_time_1(2,n_rep_each_adv_ratio) = Time_M1_to_M2_conflux;
                reaction_time_2(2,n_rep_each_adv_ratio) = time_stamp;
            end
            
        end % conflux reaction time (for n_blc)

    end% for of n_rep

    avg_reaction_time_1_cube = mean(reaction_time_1(1,:));
    avg_reaction_time_1_conflux = mean(reaction_time_1(2,:));
    
    avg_reaction_time_2_cube = mean(reaction_time_2(1,:));
    avg_reaction_time_2_conflux = mean(reaction_time_2(2,:));
    
    avg_reaction_1(1,n_adv_ratio) = avg_reaction_time_1_cube;
    avg_reaction_1(2,n_adv_ratio) = avg_reaction_time_1_conflux;
    avg_reaction_2(1,n_adv_ratio) = avg_reaction_time_2_cube;
    avg_reaction_2(2,n_adv_ratio) = avg_reaction_time_2_conflux;
end












