% decide which chain a honeset/adversarial node will choose
% under liveness attack 
% but when honest nodes did NOT notice liveness attack yet
  
function [chain_A, chain_B, weight_A, weight_B] = choose_chain_LA(new_blc, chain_A, chain_B, weight_A, weight_B)

if new_blc.type == -1 % adv %(1:honest , -1: adv) 
    % adv attaches to 2nd longest chain
    if weight_A < weight_B
        if isempty(chain_A)
            chain_A = [new_blc];
        else
            chain_A = [chain_A, new_blc];
        end
        weight_A = weight_A + new_blc.weight;
    elseif weight_A > weight_B
        if isempty(chain_B)
            chain_B = [new_blc];
        else
            chain_B = [chain_B, new_blc];
        end
        weight_B = weight_B + new_blc.weight;
    else % (weight_A == weight_B)
        if randn() > 0
            chain_B = [chain_B, new_blc];
            weight_B = weight_B + new_blc.weight;
        else
            chain_A = [chain_A, new_blc];
            weight_A = weight_A + new_blc.weight;
        end
    end
else % honest %(1:honest , -1: adv) 
    % hon attaches to 1st longest chain
    if rand() > 0.5
        if isempty(chain_A)
            chain_A = [new_blc];
        else
            chain_A = [chain_A, new_blc];
        end
        weight_A = weight_A + new_blc.weight;
    else
        if isempty(chain_B)
            chain_B = [new_blc];
        else
            chain_B = [chain_B, new_blc];
        end
        weight_B = weight_B + new_blc.weight;
    end
       
end