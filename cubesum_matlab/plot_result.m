%% plot the weight of chains

figure(1)
plot(timestamp_A_conflux, weight_A_conflux, 'b-','LineWidth',2); hold on;
plot(timestamp_B_conflux, weight_B_conflux, 'r-', 'LineWidth',2); hold on;
xline(Time_M1_to_M2_conflux,'g--','LineWidth',2); hold on;
xline(Time_M2_to_M1_conflux,'m--','LineWidth',2); hold off;
legend('Subtree 1', 'Subtree 1', 'Adaptive Weight Triggered', 'Adaptive Weight End')
title('Conflux under Liveness Attack');
xlabel('Time (sec)');
ylabel('Total weight of chain');
xlim([1,3000]);
ylim([0,6000]);

figure(2)
plot(timestamp_A_cubesum, weight_A_cubesum, 'b','LineWidth',2); hold on;
plot(timestamp_B_cubesum, weight_B_cubesum, 'r-','LineWidth',2); hold on;
xline(Time_M1_to_M2_cube,'g--','LineWidth',2); hold on;
xline(Time_M2_to_M1_cube,'m--','LineWidth',2); hold off;
legend('Subtree 1', 'Subtree 1', 'Adaptive Weight Triggered', 'Adaptive Weight End')
title('CubeSum under Liveness Attack');
xlabel('Time (sec)');
ylabel('Total weight of chain');
xlim([1,3000]);
ylim([0,6000]);
