rng(0);

% Task 1
Pt = 10;
f = 2.4 * (10^9);
lambda = 3*(10^8) / f;
Gt = 1;
Gr = 1;
d = 50;   % change distance here
Pr = Pt + 10*(log(Gt)/log(10)) +  10*(log(Gr)/log(10)) + 20*( log(lambda/(4*pi*d)) / log(10));


% Task 2
Pr_watt = 10 ^ (Pr/10 -3);
a = sqrt(Pr_watt/2) * randn(1);
b = sqrt(Pr_watt/2) * randn(1);
h = a + 1i * b;
hp = sqrt(a^2 + b^2);


% Task 3
message = randi([0, 1], 1, 300000);
% disp(message(1:10));

% BPSK
bpsk_sym = zeros(1, 300000);
for i = 1:300000
    if message(i) == 0
       bpsk_sym(i) = 1;
    elseif message(i) == 1
       bpsk_sym(i) = -1;
    end
end

% QPSK
qpsk_sym = zeros(1, 150000);
con_qpsk = [complex(1/sqrt(2), 1/sqrt(2)),
            complex(1/sqrt(2), -1/sqrt(2)),
            complex(-1/sqrt(2), 1/sqrt(2)),
            complex(-1/sqrt(2), -1/sqrt(2))];
for i = 1:150000
    mes = message(2*i-1:2*i);
    if mes(1) == 0 && mes(2) == 0
       qpsk_sym(i) = complex(1/sqrt(2), 1/sqrt(2));
    elseif mes(2) == 0 && mes(1) == 1
       qpsk_sym(i) = complex(1/sqrt(2), -1/sqrt(2));
    elseif mes(2) == 1 && mes(1) == 0
       qpsk_sym(i) = complex(-1/sqrt(2), 1/sqrt(2));
    elseif mes(1) == 1 && mes(2) == 1
       qpsk_sym(i) = complex(-1/sqrt(2), -1/sqrt(2));
    end
end   

% 16QAM
qam16_sym = zeros(1, 75000);
con_16qam = [ 3/sqrt(10) + 3/sqrt(10)*1i,
              3/sqrt(10) + 1/sqrt(10)*1i,
              3/sqrt(10) + -1/sqrt(10)*1i,
              3/sqrt(10) + -3/sqrt(10)*1i,
              1/sqrt(10) + 3/sqrt(10)*1i,
              1/sqrt(10) + 1/sqrt(10)*1i,
              1/sqrt(10) + -1/sqrt(10)*1i,
              1/sqrt(10) + -3/sqrt(10)*1i,
              -1/sqrt(10) + 3/sqrt(10)*1i,
              -1/sqrt(10) + 1/sqrt(10)*1i,
              -1/sqrt(10) + -1/sqrt(10)*1i,
              -1/sqrt(10) + -3/sqrt(10)*1i,
              -3/sqrt(10) + 3/sqrt(10)*1i,
              -3/sqrt(10) + 1/sqrt(10)*1i,
              -3/sqrt(10) + -1/sqrt(10)*1i,
              -3/sqrt(10) + -3/sqrt(10)*1i];
for i = 1:75000
    mes = message(4*i-3:4*i);
    index = mes(1) + mes(2)*2 + mes(3)*4 + mes(4)*8 +1;
    qam16_sym(i) = con_16qam(index);
end 
% disp(qam16_sym(1:3))

% 64QAM
qam64_sym = zeros(1, 50000);
con_64qam = [ 7/sqrt(42) + 7/sqrt(42)*1i,
              7/sqrt(42) + 5/sqrt(42)*1i,
              7/sqrt(42) + 3/sqrt(42)*1i,
              7/sqrt(42) + 1/sqrt(42)*1i,
              7/sqrt(42) + -1/sqrt(42)*1i,
              7/sqrt(42) + -3/sqrt(42)*1i,
              7/sqrt(42) + -5/sqrt(42)*1i,
              7/sqrt(42) + -7/sqrt(42)*1i,
              5/sqrt(42) + 7/sqrt(42)*1i,
              5/sqrt(42) + 5/sqrt(42)*1i,
              5/sqrt(42) + 3/sqrt(42)*1i,
              5/sqrt(42) + 1/sqrt(42)*1i,
              5/sqrt(42) + -1/sqrt(42)*1i,
              5/sqrt(42) + -3/sqrt(42)*1i,
              5/sqrt(42) + -5/sqrt(42)*1i,
              5/sqrt(42) + -7/sqrt(42)*1i,
              3/sqrt(42) + 7/sqrt(42)*1i,
              3/sqrt(42) + 5/sqrt(42)*1i,
              3/sqrt(42) + 3/sqrt(42)*1i,
              3/sqrt(42) + 1/sqrt(42)*1i,
              3/sqrt(42) + -1/sqrt(42)*1i,
              3/sqrt(42) + -3/sqrt(42)*1i,
              3/sqrt(42) + -5/sqrt(42)*1i,
              3/sqrt(42) + -7/sqrt(42)*1i,
              1/sqrt(42) + 7/sqrt(42)*1i,
              1/sqrt(42) + 5/sqrt(42)*1i,
              1/sqrt(42) + 3/sqrt(42)*1i,
              1/sqrt(42) + 1/sqrt(42)*1i,
              1/sqrt(42) + -1/sqrt(42)*1i,
              1/sqrt(42) + -3/sqrt(42)*1i,
              1/sqrt(42) + -5/sqrt(42)*1i,
              1/sqrt(42) + -7/sqrt(42)*1i,
              -1/sqrt(42) + 7/sqrt(42)*1i,
              -1/sqrt(42) + 5/sqrt(42)*1i,
              -1/sqrt(42) + 3/sqrt(42)*1i,
              -1/sqrt(42) + 1/sqrt(42)*1i,
              -1/sqrt(42) + -1/sqrt(42)*1i,
              -1/sqrt(42) + -3/sqrt(42)*1i,
              -1/sqrt(42) + -5/sqrt(42)*1i,
              -1/sqrt(42) + -7/sqrt(42)*1i,
              -3/sqrt(42) + 7/sqrt(42)*1i,
              -3/sqrt(42) + 5/sqrt(42)*1i,
              -3/sqrt(42) + 3/sqrt(42)*1i,
              -3/sqrt(42) + 1/sqrt(42)*1i,
              -3/sqrt(42) + -1/sqrt(42)*1i,
              -3/sqrt(42) + -3/sqrt(42)*1i,
              -3/sqrt(42) + -5/sqrt(42)*1i,
              -3/sqrt(42) + -7/sqrt(42)*1i,
              -5/sqrt(42) + 7/sqrt(42)*1i,
              -5/sqrt(42) + 5/sqrt(42)*1i,
              -5/sqrt(42) + 3/sqrt(42)*1i,
              -5/sqrt(42) + 1/sqrt(42)*1i,
              -5/sqrt(42) + -1/sqrt(42)*1i,
              -5/sqrt(42) + -3/sqrt(42)*1i,
              -5/sqrt(42) + -5/sqrt(42)*1i,
              -5/sqrt(42) + -7/sqrt(42)*1i,
              -7/sqrt(42) + 7/sqrt(42)*1i,
              -7/sqrt(42) + 5/sqrt(42)*1i,
              -7/sqrt(42) + 3/sqrt(42)*1i,
              -7/sqrt(42) + 1/sqrt(42)*1i,
              -7/sqrt(42) + -1/sqrt(42)*1i,
              -7/sqrt(42) + -3/sqrt(42)*1i,
              -7/sqrt(42) + -5/sqrt(42)*1i,
              -7/sqrt(42) + -7/sqrt(42)*1i
              ];
for i = 1:50000
    mes = message(6*i-5:6*i);
    index = mes(1) + mes(2)*2 + mes(3)*4 + mes(4)*8 + mes(5)*16 + mes(6)*32 + 1;
    qam64_sym(i) = con_64qam(index);
end   


% Task 4
n0_watt = 10 ^ ((-90/10) -3);
k = 300000;
real_noise = sqrt(n0_watt/2) * randn(1,k);
image_noise = sqrt(n0_watt/2) * randn(1,k);
noise = real_noise + 1i* image_noise;
np = mean(abs(noise).^2);
bpsk_t = bpsk_sym .* h  + noise(1:300000);
qpsk_t = qpsk_sym .* h + noise(1:150000);
qam16_t = qam16_sym .* h + noise(1:75000);
qam64_t = qam64_sym .* h + noise(1:50000);
% disp('t');
% disp(bpsk_t(1:5));
% disp('sym')
% disp(bpsk_sym(1:5).*h)
% disp('noise')
% disp(noise(1:5))


% Task 5
bpsk_r = bpsk_t ./ h;
qpsk_r = qpsk_t ./ h;
qam16_r = qam16_t ./ h;
qam64_r = qam64_t ./ h;

% BPSK
bpsk_demo = zeros(1,300000);
for i=1:300000    
    if real(bpsk_r(i))>0
        bpsk_demo(i) = 0;
    else
        bpsk_demo(i) = 1;
    end
end
% disp(bpsk_demo(1:10));

% QPSK
qpsk_demo = zeros(1,300000);
for i=1:150000
    mes = qpsk_r(i);
    if real(mes)>0
        if imag(mes)<0
            qpsk_demo(2*i-1) = 1;
        end
    else
        qpsk_demo(2*i) = 1;
        if imag(mes)<0
            qpsk_demo(2*i-1) = 1;
        end
    end
end
% disp(qpsk_demo(1:10))

% 16QAM
qam16_demo = zeros(1,300000);
for i=1:75000
    mes = qam16_r(i);
    eu_dis = abs(mes - con_16qam);
    [~, index] = min(eu_dis);
    index = index-1;
    bit = zeros(1,4);
    bit(1) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(2) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(3) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(4) = index;
    qam16_demo(4*i-3:4*i) = bit;
end
% disp(qam16_demo(1:10));

% 64QAM
qam64_demo = zeros(1,300000);
for i=1:50000
    mes = qam64_r(i);
    eu_dis = abs(mes - con_64qam);
    [~, index] = min(eu_dis);
    index = index-1;
    bit = zeros(1,6);
    bit(1) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(2) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(3) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(4) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(5) = rem(index,2);
    index = (index-rem(index,2))/2;
    bit(6) = index;
    qam64_demo(6*i-5:6*i) = bit;
end
% disp(qam64_demo(1:10));


% Task 6
% empirical noise
bpsk_n = (bpsk_r - bpsk_sym);
qpsk_n = (qpsk_r - qpsk_sym);
qam16_n = (qam16_r - qam16_sym);
qam64_n = (qam64_r - qam64_sym);
% average empirical noise power
bpsk_np = sqrt(real(bpsk_n).^2 + imag(bpsk_n).^2);
bpsk_avgnp = mean(bpsk_np);
bpsk_avgnp_dbm = 10 * (log(bpsk_avgnp*1000)/log(10));
qpsk_np = sqrt(real(qpsk_n).^2 + imag(qpsk_n).^2);
qpsk_avgnp = mean(qpsk_np);
qpsk_avgnp_dbm = 10 * (log(qpsk_avgnp*1000)/log(10));
qam16_np = sqrt(real(qam16_n).^2 + imag(qam16_n).^2);
qam16_avgnp = mean(qam16_np);
qam16_avgnp_dbm = 10 * (log(qam16_avgnp*1000)/log(10));
qam64_np = sqrt(real(qam64_n).^2 + imag(qam64_n).^2);
qam64_avgnp = mean(qam64_np);
qam64_avgnp_dbm = 10 * (log(qam64_avgnp*1000)/log(10));
% SNR
theo_snr_db = Pr - (-90);
% BPSK
xp = sqrt(real(bpsk_sym).^2 + imag(bpsk_sym).^2);
avgp = mean(xp);
avgp_dbm = 10 * (log(1000*avgp) / log(10));
bpsk_snr = avgp / bpsk_avgnp;
bpsk_snr_db = avgp_dbm - bpsk_avgnp_dbm;
% QPSK
xp = sqrt(real(qpsk_sym).^2 + imag(qpsk_sym).^2);
avgp = mean(xp);
avgp_dbm = 10 * (log(1000*avgp) / log(10));
qpsk_snr = avgp / qpsk_avgnp;
qpsk_snr_db = avgp_dbm - qpsk_avgnp_dbm;
% 16QAM
xp = sqrt(real(qam16_sym).^2 + imag(qam16_sym).^2);
avgp = mean(xp);
avgp_dbm = 10 * (log(1000*avgp) / log(10));
qam16_snr = avgp / qam16_avgnp;
qam16_snr_db = avgp_dbm - qam16_avgnp_dbm;
% 64QAM
xp = sqrt(real(qam64_sym).^2 + imag(qam64_sym).^2);
avgp = mean(xp);
avgp_dbm = 10 * (log(1000*avgp) / log(10));
qam64_snr = avgp / qam64_avgnp;
qam64_snr_db = avgp_dbm - qam64_avgnp_dbm;


% Task 7
% BER
bpsk_e = 0;
qpsk_e = 0;
qam16_e = 0;
qam64_e = 0;
for i = 1:300000
    if message(i) ~= bpsk_demo(i)
        bpsk_e = bpsk_e + 1;
    end
    if message(i) ~= qpsk_demo(i)
        qpsk_e = qpsk_e + 1;
    end
    if message(i) ~= qam16_demo(i)
        qam16_e = qam16_e + 1;
    end
    if message(i) ~= qam64_demo(i)
        qam64_e = qam64_e + 1;
    end
end
bpsk_ber = bpsk_e / 300000;
qpsk_ber = qpsk_e / 300000;
qam16_ber = qam16_e / 300000;
qam64_ber = qam64_e / 300000;
% theoretical PDR
packet_bit = 4000;
bpsk_theo_pdr = (1-bpsk_ber)^packet_bit;
qpsk_theo_pdr = (1-qpsk_ber)^packet_bit;
qam16_theo_pdr = (1-qam16_ber)^packet_bit;
qam64_theo_pdr = (1-qam64_ber)^packet_bit;
% empirical PDR
bpsk_er_packet = 0;
qpsk_er_packet = 0;
qam16_er_packet = 0;
qam64_er_packet = 0;
for j = 1:75
    for i = (packet_bit*j-(packet_bit-1)):(packet_bit*j)
        if message(i) ~= bpsk_demo(i)
           bpsk_er_packet = bpsk_er_packet+1;
           break;
        end
    end
end
bpsk_pdr = 1 - (bpsk_er_packet / (300000/packet_bit));
for j = 1:75
    for i = (packet_bit*j-(packet_bit-1)):(packet_bit*j)
        if message(i) ~= qpsk_demo(i)
           qpsk_er_packet = qpsk_er_packet+1;
           break;
        end
    end
end
qpsk_pdr = 1 - (qpsk_er_packet / (300000/packet_bit));
for j = 1:75
    for i = (packet_bit*j-(packet_bit-1)):(packet_bit*j)
        if message(i) ~= qam16_demo(i)
           qam16_er_packet = qam16_er_packet+1;
           break;
        end
    end
end
qam16_pdr = 1 - (qam16_er_packet / (300000/packet_bit));
for j = 1:75
    for i = (packet_bit*j-(packet_bit-1)):(packet_bit*j)
        if message(i) ~= qam64_demo(i)
           qam64_er_packet = qam64_er_packet+1;
           break;
        end
    end
end
qam64_pdr = 1 - (qam64_er_packet / (300000/packet_bit));
% theoretical throughput
time = 3.2*10^(-6);   % sample duration
bpsk_theo_tput = bpsk_theo_pdr / time;
qpsk_theo_tput = qpsk_theo_pdr * 2 / time;
qam16_theo_tput = qam16_theo_pdr * 4 / time;
qam64_theo_tput = qam64_theo_pdr * 6 / time;
% empirical throughput
bpsk_tput = bpsk_pdr  / time;
qpsk_tput = qpsk_pdr * 2 / time;
qam16_tput = qam16_pdr * 4 / time;
qam64_tput = qam64_pdr * 6 / time;
% disp
disp(['distance:',num2str(d)]);
disp(['Prx (Watt): ',num2str(Pr_watt)]);
disp(['Prx (dbm): ',num2str(Pr)]);
disp(['theoretical noise power:',num2str(10^(-12)),'(Watt) = -90(dbm)'])
disp(['theoretical SNR:', num2str(10^(theo_snr_db/10))])
disp(['theoretical SNR(dB):', num2str(theo_snr_db)])


disp('BPSK:');
disp(['empirical average noise power:',num2str(bpsk_avgnp), '(Watt) = ', num2str(bpsk_avgnp_dbm), '(dbm)'])
disp(['empirical SNR:',num2str(bpsk_snr)]);
disp(['empirical SNR(dB):',num2str(bpsk_snr_db)]);
disp(['empirical BER:',num2str(bpsk_ber)]);
disp(['empirical throughput:',num2str(bpsk_tput),' (bit/s)']);
disp(['theoretical throughput:',num2str(bpsk_theo_tput),' (bit/s)']);
disp('QPSK:');
disp(['empirical average noise power:',num2str(qpsk_avgnp), '(Watt) = ', num2str(qpsk_avgnp_dbm), '(dbm)'])
disp(['empirical SNR:',num2str(qpsk_snr)]);
disp(['empirical SNR(dB):',num2str(qpsk_snr_db)]);
disp(['empirical BER:',num2str(qpsk_ber)]);
disp(['empirical throughput:',num2str(qpsk_tput),' (bit/s)']);
disp(['theoretical throughput:',num2str(qpsk_theo_tput),' (bit/s)']);
disp('16QAM:');
disp(['empirical average noise power:',num2str(qam16_avgnp), '(Watt) = ', num2str(qam16_avgnp_dbm), '(dbm)'])
disp(['empirical SNR:',num2str(qam16_snr)]);
disp(['empirical SNR(dB):',num2str(qam16_snr_db)]);
disp(['empirical BER:',num2str(qam16_ber)]);
disp(['empirical throughput:',num2str(qam16_tput),' (bit/s)']);
disp(['theoretical throughput:',num2str(qam16_theo_tput),' (bit/s)']);
disp('64QAM:');
disp(['empirical average noise power:',num2str(qam64_avgnp), '(Watt) = ', num2str(qam64_avgnp_dbm), '(dbm)'])
disp(['empirical SNR:',num2str(qam16_snr)]);
disp(['empirical SNR(dB):',num2str(qam64_snr_db)]);
disp(['empirical BER:',num2str(qam64_ber)]);
disp(['empirical throughput:',num2str(qam64_tput),' (bit/s)']);
disp(['theoretical throughput:',num2str(qam64_theo_tput),' (bit/s)']);
disp(' ')
em_tput = [bpsk_tput, qpsk_tput, qam16_tput, qam64_tput];
[optimal_tput, index] = max(em_tput);
if index == 1
    disp('optimal modulation scheme: BPSK')
elseif index ==2
    disp('optimal modulation scheme: QPSK')
elseif index == 3
    disp('optimal modulation scheme: 16QAM')
else
    disp('optimal modulation scheme: 64QAM')
end

% Task 8
% plot_sig = bpsk_r;
% plot_sig = qpsk_r;
% plot_sig = qam16_r;
% plot_sig = qam64_r;
% scatter(real(plot_sig), imag(plot_sig), '.');
% xlabel('I');
% ylabel('Q');
% xlim([-1.5, 1.5]);
% ylim([-1.5, 1.5]);
% title('BPSK Costellation Diagram')
% title('QPSK Costellation Diagram')
% title('16QAM Costellation Diagram')
% title('64QAM Costellation Diagram')
% grid on;


% x = 50:50:600;
% y_bpsk = [312500, 312500, 312500, 312500, 312500, 312500, 312500, 308333.3333, 287500, 229166.6667, 100000, 8333.3333];
% y_qpsk = [625000, 625000, 625000, 625000, 625000, 608333.3333, 475000, 66666.6667, 0, 0, 0, 0];
% y_16qam = [1250000, 1250000, 1216666.6667, 33333.3333, 0, 0, 0, 0, 0, 0, 0, 0];
% y_64qam = [1875000, 50000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% plot(x, y_bpsk, x, y_qpsk, x, y_16qam, x, y_64qam);
% xlabel('Distance (m)');
% ylabel('Throughput (bit/s)');
% title('Throughput vs Distance');
% legend('BPSK', 'QPSK', '16QAM', '64QAM');
% grid on;




% Question
load('SNR_BER.mat')
bpsk_snrdb_round = round(bpsk_snr_db);
qpsk_snrdb_round = round(qpsk_snr_db);
qam16_snrdb_round = round(qam16_snr_db);
qam64_snrdb_round = round(qam64_snr_db);

bpsk_theo_ber = SNR_BER(1,bpsk_snrdb_round);
qpsk_theo_ber = SNR_BER(2,qpsk_snrdb_round);
qam16_theo_ber = SNR_BER(3,qam16_snrdb_round);
qam64_theo_ber = SNR_BER(4,qam64_snrdb_round);

packet_bit = 100;
bpsk_table_pdr = (1-bpsk_theo_ber)^packet_bit;
qpsk_table_pdr = (1-qpsk_theo_ber)^packet_bit;
qam16_table_pdr = (1-qam16_theo_ber)^packet_bit;
qam64_table_pdr = (1-qam64_theo_ber)^packet_bit;

bpsk_table_tput = bpsk_table_pdr / time;
qpsk_table_tput = qpsk_table_pdr * 2 / time;
qam16_table_tput = qam16_table_pdr * 4 / time;
qam64_table_tput = qam64_table_pdr * 6 / time;

theo_tput = [bpsk_table_tput, qpsk_table_tput, qam16_table_tput, qam64_table_tput];
[optimal_tput, index] = max(theo_tput);
if index == 1
    disp('theoretical optimal modulation scheme: BPSK')
elseif index ==2
    disp('theoretical optimal modulation scheme: QPSK')
elseif index == 3
    disp('theoretical optimal modulation scheme: 16QAM')
else
    disp('theoretical optimal modulation scheme: 64QAM')
end

