% decide which chain a honeset/adversarial node will choose
% under liveness attack 
  
function [flag_cube, Delta_out, B_out] = cube_sum_check(chain_A, chain_B, B_in, Delta_threshold, B_threshold, weight_A, weight_B)

if weight_A > weight_B
    heavier_weight = weight_A;
    lighter_weight = weight_B;
else
    heavier_weight = weight_B;
    lighter_weight = weight_A;
end

Delta_out = abs(weight_A - weight_B);

if heavier_weight-lighter_weight <= Delta_threshold
    B_out = B_in +1; % increment B by 1
else
    B_out = 0; % reset B
end

if B_out >= B_threshold
   flag_cube = 1; %1: mode 2, 0 : mode 1
else
    flag_cube = 0;
end
