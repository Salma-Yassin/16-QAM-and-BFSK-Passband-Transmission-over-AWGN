%% Setting parameters 
clear 
clc 

Tb = 5; %% bit duration in secs 
T = 4 * Tb;  %%quadbit duration 

E0 = 1; %% enerdy of the signal with the smallest amplitude 
M = 16;
Eav = 2/3 * (M-1) * E0;

N_bit = 500; %%number of samples per bit 
N_quadbit = 4 * N_bit;

t_bit = linspace(0, Tb ,N_bit); %% time base for each bit  
t_quadbit = linspace(0, T , N_quadbit);

msg_l = 36 ; %% number of bits sent which has to be divisible by 4---> 300 point 
t_signal = linspace(0, msg_l*Tb, msg_l*N_bit);%% total duration of the messag 

fc = 2/Tb; %% frequency of the carrier 
N0 = 5;

%% Message source ------> generates a rondom stream of 0s and 1s

message = randi([0 1],1,msg_l);
figure
stem(message)
title('Randomly generated Digital message')
ylabel('Bit value')
xlabel('Bit Index')
ylim([-0.5 1.5])


%% Signal transimission encoder --------> polar non return tozero: (1)->1,(0)->-1

encodedMessage = [];

for i = 1:4:length(message)
    d = bi2de(flip(message(i:i+3)));
    encodedMessage = [encodedMessage d];
end 

a = [-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,1,1,1,1]*sqrt(E0);
b = [-3,-1,3,1,-3,-1,3,1,-3,-1,3,1,-3,-1,3,1]*sqrt(E0);

%% 16-QAM modulator --->  

carrier_i = sqrt(2/T) * cos(2 * pi * fc * t_quadbit); %% normalize basis function 
carrier_Q = sqrt(2/T) * sin(2 * pi * fc * t_quadbit); %% normalize basis function 

modulatedSignal = [];

for i = 1:1:length(encodedMessage)
    index = encodedMessage(i) + 1;
    signal = a(index)* carrier_i + b(index)* carrier_Q;
    modulatedSignal = [modulatedSignal signal];
end 


plot(t_signal,modulatedSignal)
title('Modulated 16-QAM signal')
ylabel('Amplitude of the modulated signal')
xlabel('Time in secs')
grid on
legend('16-QAM signal')


%% Constellation of transsmitted 16-QAM

basis_func_i = sqrt(2/T) * cos(2 * pi * fc * t_quadbit);
basis_func_Q = sqrt(2/T) * sin(2 * pi * fc * t_quadbit);

si1_vector=[];
si2_vector=[];

for i = 1:N_quadbit:length(modulatedSignal)
    vec = modulatedSignal(i: i + N_quadbit - 1);
    vec = vec.*basis_func_i;
    intg = trapz(t_quadbit,vec); %% seperation is tb 
    si1_vector = [si1_vector intg];
end  

for i = 1:N_quadbit:length(modulatedSignal)
    vec = modulatedSignal(i: i + N_quadbit - 1);
    vec = vec.*basis_func_Q;
    intg = trapz(t_quadbit,vec); %% seperation is tb 
    si2_vector = [si2_vector intg];
end 

si_vector=[si1_vector ; si2_vector];
scatterplot(transpose(si_vector));
title('Constellation of the Transimitted 16-QAM')
ylabel('Q Phase basis')
xlabel('I phase basis')
hold on
x = [-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,1,1,1,1]*sqrt(E0);
y = [-3,-1,3,1,-3,-1,3,1,-3,-1,3,1,-3,-1,3,1]*sqrt(E0);
labels = ['0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'];
grid on
for i = 1:1:length(x)
    %text(x(i),y(i),labels(i),'VerticalAlignment','bottom','HorizontalAlignment','right');
    index = [0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60];
    index = index + 1;
    text(x(i),y(i),labels(index(i) : index(i) + 3),'VerticalAlignment','bottom','HorizontalAlignment','right');
end 


%% Adding whie additive Gussian noise from a normal distribution 

recievedSignal = modulatedSignal + unifrnd(0,N0/2,1,length(modulatedSignal));

figure
plot(t_signal,recievedSignal)
title(' Recieved Modulated 16-QAM signal with AWGN')
ylabel('Amplitude of Recieved modulated signal')
xlabel('Time in secs')
grid on
xlim([0 50])
legend('16-QAM signal')

%% Constellation of recieved 16-QAM

basis_func_i = sqrt(2/T) * cos(2 * pi * fc * t_quadbit);
basis_func_Q = sqrt(2/T) * sin(2 * pi * fc * t_quadbit);

xi1_vector=[];
xi2_vector=[];

for i = 1 : N_quadbit : length(recievedSignal)
    vec = recievedSignal(i: i + N_quadbit - 1);
    vec = vec.*basis_func_i;
    intg = trapz(t_quadbit ,vec); %% seperation is tb 
    xi1_vector = [xi1_vector intg];
end  

for i = 1 : N_quadbit : length(recievedSignal)
    vec = recievedSignal(i: i + N_quadbit - 1);
    vec = vec.*basis_func_Q;
    intg = trapz(t_quadbit,vec); %% seperation is tb 
    xi2_vector = [xi2_vector intg];
end 

xi_vector=[xi1_vector ; xi2_vector];
scatterplot(transpose(xi_vector))
title('Constellation of the Received 16-QAM')
ylabel('Q Phase basis')
xlabel('I phase basis')
grid on

hold on
x = [-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,1,1,1,1]*sqrt(E0);
y = [-3,-1,3,1,-3,-1,3,1,-3,-1,3,1,-3,-1,3,1]*sqrt(E0);
labels = ['0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'];
grid on
for i = 1:1:length(x)
    %text(x(i),y(i),labels(i),'VerticalAlignment','bottom','HorizontalAlignment','right');
    index = [0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60];
    index = index + 1;
    text(x(i),y(i),labels(index(i) : index(i) + 3),'VerticalAlignment','bottom','HorizontalAlignment','right');
end 

%% Signal transimission decoder -----> using ML rule

% calculate the minimum distance of the given point
point_x = [-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,1,1,1,1]*sqrt(E0);
point_y = [-3,-1,3,1,-3,-1,3,1,-3,-1,3,1,-3,-1,3,1]*sqrt(E0);
codes = [0,0,0,0 ; 0,0,0,1 ; 0,0,1,0 ; 0,0,1,1 ; 0,1,0,0
        ;0,1,0,1 ; 0,1,1,0 ; 0,1,1,1 ; 1,0,0,0 ; 1,0,0,1 
        ; 1,0,1,0 ; 1,0,1,1 ; 1,1,0,0 ; 1,1,0,1 ; 1,1,1,0 
        ; 1,1,1,1];

rec_signal_Decoded = [];

for i = 1:1:length(xi_vector)
    dmin_point = inf;
    index_point = 0;
    for j = 1:1:16 % iterate over the 16 point to get the min distance 
        d = sqrt((xi_vector(1,i)- point_x(j) )^2 + (xi_vector(2,i)- point_y(j))^2);
        if d < dmin_point
            dmin_point = d;
            index_point = j;
        end 
    end 
    rec_signal_Decoded = [rec_signal_Decoded codes(index_point,:)];
end 

%% check

message ;
rec_signal_Decoded;
ans = message - rec_signal_Decoded;

%% Psd of the baseband 16-QAM 

base_band = [];

for i = 1:1:length(encodedMessage)
    index = encodedMessage(i) + 1;
    seg = a(index)*ones(1,N_quadbit) + b(index)* ones(1,N_quadbit);
    base_band = [base_band seg];
end 

ts = Tb/N_bit;
[psd,f]= periodogram(base_band,[],[],1/ts*50);

f = f/Tb;

psd = 10*log10(psd/(2*E0));
%psd = filter([ 0.5 0.5 0.5 0.5 0.5],1,psd);
% figure
% %plot(f,psd)
% xlim([0 5])
% title('PSD of the Transimitted 16-QAM')
% ylabel('Spectral power Density in dB ')
% xlabel('Normalized Frequency')
% grid on
% legend('16-QAM PSD')



