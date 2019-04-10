%% initialization

format short; format compact;
clear all; clc; close all;

%% load data
%%%%%%% Load File %%%%%%%%
filename    = 'fm_5kbps_d06_batteryfree_wired.lvm';
rx          = lvmread(filename);
rx          = rx(:,2);

%%%%%%% set 1's and 0's %%%%%%%%
th = max(rx)/2;
for i = 1:length(rx)
    if rx(i) >= th
        rx(i) = 1;
    else
        rx(i) = 0;
    end
end

%%%%%%%% Load Message %%%%%%%%
msg_s = load('bf_data.mat');
msg = transpose(msg_s.data);


%% reset
%%%%%%%% clear vars for testing purposes %%%%%%%%
clearvars -except filename rx th msg_s msg
close all

%% data cleaning

%%%%%%%% remove blips of 1's and 0's %%%%%%%%
clean_count = [];

for i = 2:length(rx)
    % if 0,1,0
    if rx(i-1) == 0 & rx(i) == 1 & rx(i+1) == 0
        rx_c(i) = 0;
        clean_count = [clean_count,0];
    % if 1,0,1
    elseif rx(i-1) == 1 & rx(i) == 0 & rx(i+1) == 1
        rx_c(i) = 1;
        clean_count = [clean_count,1];
    else
        rx_c(i) = rx(i);
    end
end
rx_c(1) = 1;
fprintf('number of bits cleaned: %d\n',length(clean_count))


%% data processing

%%%%%%%% Set Parameters %%%%%%%%
spb = 384;      % samples per bit
rx_d = [];      % discretized rx
k_check = [];   % bit length matrix, debug
j = 1;          % signal indexer
k=1;            % bit counter
m=0.75245;      % slope of correction
b=95.057;       % y-int of correction

%%%%%%%% discretize signal %%%%%%%%
for i = 2:length(rx_c)
    % if bits are the same
    if rx_c(i) - rx_c(i-1) == 0 
        k=k+1;
    % if 0 -> 1
    elseif rx_c(i) - rx_c(i-1) == 1
        n = round((m*k+b)/spb);         % correction function
%         k_check(j:j+n,1) = i;           % debug
%         k_check(j:j+n,2) = k;           % debug
%         k_check(j:j+n,3) = n;           % debug
        rx_d(j:j+n) = 0;
        j = j + n;
        k = 1;
    % if 1 -> 0
    elseif rx_c(i) - rx_c(i-1) == -1
        n = round((m*k+b)/spb);         % correction function
%         k_check(j:j+n,1) = i;           % debug
%         k_check(j:j+n,2) = k;           % debug
%         k_check(j:j+n,3) = n;           % debug
        rx_d(j:j+n) = 1;
        j = j + n;
        k = 1;
    end
end

rx_n = rx_d;


%%%%%%%% debugging matrices %%%%%%%%

% thanks k_checku and k_tab, they are the MVMs
k_checku = unique(k_check,'rows','stable');
k_tabd = tabulate(k_check(:,2));
k_tab = k_tabd(k_tabd(:,2) > 0,:);
k_tab(:,4) = k_tab(:,1)/384;
k_tab(:,5) = k_tab(:,1)/128;

%fprintf('\nsamples in rx_c: %f\n',length(rx_c)/spb)
%fprintf('samples in rx_d: %d\n',length(rx_d))


%%%%%%%% debugging plots %%%%%%%%
% figure
% 
% plt_st = 75;
% plt_end = 110;
% subplot(2,1,1)
% stem(rx_c(plt_st*spb:plt_end*spb))
% title('clean rx (rx_c)')
% subplot(2,1,2)
% stem(rx_d(plt_st:plt_end))
% title('discrete rx (rx_d)')
% 
% figure
% kx = k_tab(1:44,4);
% ky = k_tab(1:44,2);
% 
% bar(kx,ky)
% set(gca,'yscale','log');
% xlabel('bit length divided by spb');
% ylabel('number of occurrences');
% title('spb = 384');


%% Main

%%%%%%%% create header %%%%%%%%
header = ones(10, 1); header(1:2:10) = 0; header = transpose(header);

%%%%%%%% find data by locating header %%%%%%%%
h_start = strfind(rx_n, header);

%%%%%%%% removing footers and mismatches %%%%%%%%
check = 1;

while check ~= 0
    check = 0;
    for i=2:length(h_start)        
        if h_start(i) - h_start(i-1) < 200
            h_start(i) = 0;
            check = 1;
        end
    end
    h_start = h_start(h_start ~= 0);
    fprintf('iter\n')
end

%%%%%%%% loading data from header locations %%%%%%%%
d_start = h_start + length(header);
d_end = d_start + (length(msg)-1);

for i = 1:(length(d_start)-1)
    data_rx(i,:) = rx_n(d_start(i):d_end(i));
end

%%%%%%%% find data by direct match (cheat) %%%%%%%%
d_loc = strfind(rx_n, msg.');
% d_loc_diff = (d_loc - [0,d_loc(1:end-1)]).';            % debug

for i = 1:length(d_loc)
    data_cheat(i,:) = rx_n(d_loc(i):d_loc(i)+(length(msg)-1));
end

%%%%%%%% ber %%%%%%%%
ber_count   = 1;
ber         = [];

for i = 1:size(data_rx,1)
   ber(i) = sum(abs(data_rx(i,:)-msg))/length(msg);
end

ber
ber_cum = mean(ber ~= 0)


%% Data Plotting

%%%%%%%% Data by direct match (cheat) %%%%%%%%
figure
for i = 1:size(data_cheat,1)
    hold on;
    stem(data_cheat(i,:)); 
end
axis([0 length(data_cheat) -0.1 1.1]);

%%%%%%%% Data by finding header %%%%%%%%
figure
for i = 1:size(data_rx,1)
    hold on;
    stem(data_rx(i,:)); 
end
axis([0 length(data_cheat) -0.1 1.1]);

