%% BERs plots for 16-QAM and BFSK both theortical and through simulations

clear
clc

N0=10;
Eb=[];
x=0.0001;
Eb=[Eb x];
for i=1:1:5 %generating a vector of Eb values ranging from 0.0001 to 
    x=x*10;
    for j=1:1:9
        Eb=[Eb j*x];
    end
end

%% BER Theortical

M = 16;
E0 = Eb;
%Eav = 2/3 * (M-1) * E0;
Eav = E0;

BER_16_QAM = (1/log2(M))*(2*(1-1/sqrt(M))*erfc(sqrt((3*Eav)/(2*(M-1)*N0))));
BER_BFSK = 0.5 * erfc(sqrt(Eb/(2*N0))); % BER for binary FSK 

subplot(2,1,1)
semilogy(20*log10(Eb/N0),BER_16_QAM)
hold on 
semilogy(20*log10(Eb/N0),BER_BFSK)
title('BERs of 16 - QAM and BFSK Theortical')
xlabel('Eb/N0(dB)')
ylabel('Bit Error Rate')
legend('16 - QAM','BFSK')
grid on 

%% BER through Simulation 

BER_BFSK_sim = [];
BER_16_QAM_sim=[];

for i = 1:1:length(Eb) 
    
    b = BFSK_BER(Eb(i),N0);
    BER_BFSK_sim = [BER_BFSK_sim b];
end

Eav = Eb;
E0 = 3/2*Eav/(M-1);
for i =1:1:length(E0)% run the 16-QAM modulation for different valus of E0/N0
    b = QAM_16_BER(E0(i),N0);
    BER_16_QAM_sim=[BER_16_QAM_sim b];
end

subplot(2,1,2)
semilogy(20*log10(Eb/N0),BER_16_QAM_sim)
hold on 
semilogy(20*log10(Eb/N0),BER_BFSK_sim)

title('BERs of 16-QAM and BFSK Simulation')
xlabel('Eb/N0(dB)')
ylabel('Bit Error Rate')
legend('16 - QAM','BFSK')
grid on 


%% Simulation Errors of BPSK , QPSK, 16-QAM, BFSK

BER_BPSK_sim=[];
BER_QPSK_sim=[];
BER_BFSK_sim = [];
BER_16_QAM_sim=[];

for i =1:1:length(Eb) % run the BPSK modulation for different valus of Eb/N0
    b = BPSK_mode(N0,Eb(i));
    BER_BPSK_sim = [BER_BPSK_sim b];
end

for i =1:1:length(Eb)% run the QPSK modulation for different valus of Eb/N0
    b = QPSK_mode(N0,Eb(i));
    BER_QPSK_sim=[BER_QPSK_sim b];
end


for i = 1:1:length(Eb) % run the BFSK modulation for different valus of Eb/N0
    b = BFSK_BER(Eb(i),N0);
    BER_BFSK_sim = [BER_BFSK_sim b];
end

Eav = Eb;
E0 = 3/2*Eav/(M-1);

for i =1:1:length(Eb)% run the 16-QAM modulation for different valus of E0/N0
    b = QAM_16_BER(E0(i),N0);
    BER_16_QAM_sim=[BER_16_QAM_sim b];
end
 
%%
figure()

semilogy(20*log10(Eb/N0),BER_BPSK_sim)
hold on 
semilogy(20*log10(Eb/N0),BER_QPSK_sim)
semilogy(20*log10(Eb/N0),BER_16_QAM_sim)
semilogy(20*log10(Eb/N0),BER_BFSK_sim)

title('BERs of BPSK, QPSK, 16 - QAM, BFSK Simulation')
xlabel('Eb/N0(dB)')
ylabel('Bit Error Rate')
legend('BPSK','QPSK','16 - QAM','BFSK')
grid on 
xlim([-100 19.08])
