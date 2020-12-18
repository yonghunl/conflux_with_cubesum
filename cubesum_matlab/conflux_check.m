% decide which chain a honeset/adversarial node will choose
% under liveness attack 
  
function flag_conflux = conflux_check(chain_A, chain_B, alpha, beta, gamma, new_time_stamp, weight_A, weight_B)

heavier_weight = max(weight_A, weight_B);
lighter_weight = min(weight_A, weight_B);

% the default mode is 1, implying flag = 0 
% flag = 1 only when operating mode 2 and satisfying all the 3 conditions.
flag_conflux = 0; 

% "f(a) < alpha" <=> "subTW(a) - sibling's subTW(a) < alpha" 
if heavier_weight - lighter_weight <= alpha
    
    if (isempty(chain_A) ~= 1) && (isempty(chain_B) ~= 1)
        % "t(a) > beta"
        % t(a) = |Timerchain(b)| - |Timerchain(a.parent)| 
        % where b is the new block
        genesis_time_stamp = 0;
        if (new_time_stamp - genesis_time_stamp) >= beta
            flag_conflux = 1;

        % "g(a) > gamma"
        % g(a) is subTW(a.parent)
        % weight of chain A and B and genesis.
        elseif weight_A + weight_B + 1 > gamma 
            flag_conflux = 1;
        end
    end
end



