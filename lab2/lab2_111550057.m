close all; clear; clc;
addpath ./tasks;

% Environment Configurations
tx_node_number = 1;          % Number of Tx users
rx_node_number = 2;          % Number of Rx users
analog_antenna_number = 16;  % Number of Tx antennas of analog
digital_antenna_number = 2;  % Number of Tx antennas of digital
rx_antenna_number = 1;       % Number of Rx antennas

% Generate receivers with beam direction and distances
% [x y] = rand_rx_location_list(tx_beam_direction_idx)
rand_rx_location_list = [];
for tx_beam = 0:10:180
    tmp = [];
    for d = 50:50:500
        offset = -5 + 10 * rand();       % -5~5 degrees
        x = d * cosd(tx_beam + offset);  % Add a small random offset
        y = d * sind(tx_beam + offset);  % Add a small random offset
        tmp = [tmp x y];
    end
    rand_rx_location_list = [rand_rx_location_list; tmp];
end

% Define coordination and power for transmitter
origin = [0, 0];
tx_location = origin;
P_tx_dBm = 10;          % Transmission power of Tx (dBm)
N0_dBm = -95;           % Assume noise power is -90 dBm

% Randomly select receiver's coordination at d=200 -> (7:8) is d=200
% rand_idx = randperm(numel(0:10:180), 2);
% rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
% rx2_location = rand_rx_location_list(rand_idx(2), (7:8));

% TODO: Implement Analog Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
rx1_avgSNR_list_8 = [];
rx2_avgINR_list_8 = [];

for index = 1:2:19
    rx1_SNR_list = [];
    rx2_INR_list = [];
    for i = 1:10
        rand_idx = randperm(numel(0:10:180), 2);
        rx1_location = rand_rx_location_list(rand_idx(1), (index:index+1));
        rx2_location = rand_rx_location_list(rand_idx(2), (index:index+1));

        [rx1_SNR_dbm, rx2_INR_dbm, ~] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, 8, 0:10:180);

        rx1_SNR_list = [rx1_SNR_list rx1_SNR_dbm];
        rx2_INR_list = [rx2_INR_list rx2_INR_dbm];

        % fprintf('Task1: SNR Calculation\n');
        % fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);
    end
    rx1_avgSNR_list_8 = [rx1_avgSNR_list_8 mean(rx1_SNR_list)];
    rx2_avgINR_list_8 = [rx2_avgINR_list_8 mean(rx2_INR_list)];
end

rx1_avgSNR_list_16 = [];
rx2_avgINR_list_16 = [];

for index = 1:2:19
    rx1_SNR_list = [];
    rx2_INR_list = [];
    for i = 1:10
        rand_idx = randperm(numel(0:10:180), 2);
        rx1_location = rand_rx_location_list(rand_idx(1), (index:index+1));
        rx2_location = rand_rx_location_list(rand_idx(2), (index:index+1));

        [rx1_SNR_dbm, rx2_INR_dbm, ~] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, 16, 0:10:180);

        rx1_SNR_list = [rx1_SNR_list rx1_SNR_dbm];
        rx2_INR_list = [rx2_INR_list rx2_INR_dbm];

        % fprintf('Task1: SNR Calculation\n');
        % fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);
    end
    rx1_avgSNR_list_16 = [rx1_avgSNR_list_16 mean(rx1_SNR_list)];
    rx2_avgINR_list_16 = [rx2_avgINR_list_16 mean(rx2_INR_list)];
end

d = 50:50:500;

figure(1); 
hold on;
plot(d, rx1_avgSNR_list_8, '-r', 'DisplayName', 'R_{x1} SNR (8 antenna)'); 
plot(d, rx2_avgINR_list_8, '-g', 'DisplayName', 'R_{x2} INR (8 antenna)'); 
plot(d, rx1_avgSNR_list_16, '-b', 'DisplayName', 'R_{x1} SNR (16 antenna)'); 
plot(d, rx2_avgINR_list_16, '-k', 'DisplayName', 'R_{x2} INR (16 antenna)');
hold off;
xlabel('Distance (m)');
xlim([0 550]);
ylabel('SNR or INR (db)');
title('Average SNR/INR vs. Distance');
legend('show');
grid on;


rx1_SNR_list = [];
rx2_INR_list = [];
for i = 1:10
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));

    [rx1_SNR_dbm, rx2_INR_dbm, ~] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number, 0:10:180);

    rx1_SNR_list = [rx1_SNR_list rx1_SNR_dbm];
    rx2_INR_list = [rx2_INR_list rx2_INR_dbm];

    fprintf('Task1: SNR Calculation\n');
    fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);
end

topo = 1:10;

figure(2); 
hold on;
plot(topo, rx1_SNR_list, '-b', 'DisplayName', 'R_{x1} SNR'); 
plot(topo, rx2_INR_list, '-r', 'DisplayName', 'R_{x2} INR'); 
hold off;
xlabel('Topology no.');
ylabel('SNR or INR (db)');
title('SNR/INR of 10 Topologies');
legend('show');
grid on;


P_rx_list_19 = [];
P_rx_list_37 = [];
P_rx_list_73 = [];
for i = 1:10
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));

    [~, ~, P_rx1_dbm] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number, 0:10:180);
    P_rx_list_19 = [P_rx_list_19 P_rx1_dbm];

    [~, ~, P_rx1_dbm] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number, 0:5:180);
    P_rx_list_37 = [P_rx_list_37 P_rx1_dbm];

    [~, ~, P_rx1_dbm] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number, 0:2.5:180);
    P_rx_list_73 = [P_rx_list_73 P_rx1_dbm];

    % fprintf('Task1: SNR Calculation\n');
    % fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);
end
topo = 1:10;

figure(3); 
hold on;
plot(topo, P_rx_list_19, '-b', 'DisplayName', 'codebook size = 19', 'LineWidth',3);
plot(topo, P_rx_list_37, '-r', 'DisplayName', 'codebook size = 37', 'LineWidth',2); 
plot(topo, P_rx_list_73, '-g', 'DisplayName', 'codebook size = 73', 'LineWidth',1); 
hold off;
xlabel('Topology no.');
xlim([0 11]);
ylabel('P_{rx} (dbm)');
title('P_{rx} of 10 Topologies');
legend('show');
grid on;


% TODO: Implement Digital Beamforming functions in /tasks
% Hint: you can adjust input/output for reports

rx1_avgSNR_list = [];
rx2_avgSNR_list = [];
rx1_avgSNR_list_ori = [];
rx2_avgSNR_list_ori = [];

for index = 1:2:19
    rx1_SNR_list = [];
    rx2_SNR_list = [];
    rx1_SNR_list_ori = [];
    rx2_SNR_list_ori = [];
    for i = 1:10
        rand_idx = randperm(numel(0:10:180), 2);
        rx1_location = rand_rx_location_list(rand_idx(1), (index:index+1));
        rx2_location = rand_rx_location_list(rand_idx(2), (index:index+1));

        [rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm, ~, ~] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);

        rx1_SNR_list = [rx1_SNR_list rx1_SNR_dbm];
        rx2_SNR_list = [rx2_SNR_list rx2_SNR_dbm];
        rx1_SNR_list_ori = [rx1_SNR_list_ori ori_rx1_SNR_dbm];
        rx2_SNR_list_ori = [rx2_SNR_list_ori ori_rx2_SNR_dbm];

        fprintf('Task2: SNR Calculation\n');
        fprintf('Receiver1\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx1_SNR_dbm, ori_rx1_SNR_dbm);
        fprintf('Receiver2\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx2_SNR_dbm, ori_rx2_SNR_dbm);
    end
    rx1_avgSNR_list = [rx1_avgSNR_list mean(rx1_SNR_list)];
    rx2_avgSNR_list = [rx2_avgSNR_list mean(rx2_SNR_list)];
    rx1_avgSNR_list_ori = [rx1_avgSNR_list_ori mean(rx1_SNR_list_ori)];
    rx2_avgSNR_list_ori = [rx2_avgSNR_list_ori mean(rx2_SNR_list_ori)];
end
d = 50:50:500;

figure(4); 
hold on;
plot(d, rx1_avgSNR_list, '-r', 'DisplayName', 'R_{x1} SNR w/ ZFBF','LineWidth',2); 
plot(d, rx2_avgSNR_list, '-g', 'DisplayName', 'R_{x2} SNR w/ ZFBF','LineWidth',1); 
plot(d, rx1_avgSNR_list_ori, '-b', 'DisplayName', 'R_{x1} SNR w/o ZFBF','LineWidth',2); 
plot(d, rx2_avgSNR_list_ori, '-k', 'DisplayName', 'R_{x2} SNR w/o ZFBF','LineWidth',1);
hold off;
xlabel('Distance (m)');
xlim([0 550]);
ylabel('SNR (db)');
title('Average SNR vs. Distance');
legend('show');
grid on;

error_list = [];
H_eq_list = [];
for i = 1:10
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));

    [~, ~, ~, ~, H_eq, error] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);
    
    error_list = [error_list error];
    H_eq_list = [H_eq_list H_eq(1)];
end

figure(5); 
plot(1:10, H_eq_list, '-b','LineWidth',2); 
xlabel('Topology No.');
ylabel('H_{eq}');
title('H_{eq} of R_1 w/ ZFBF');
grid on;

figure(6); 
plot(1:10, error_list, '-b','LineWidth',2); 
xlabel('Topology No.');
ylabel('error (dbm)');
title('error of R_1 w/ ZFBF');
grid on;
