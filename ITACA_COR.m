%{ 
%% ELECTROCARDIOGRAM SIGNAL PROCESSING TEST
   _____________________________________
  |                                     |
  |Roberto Holgado Cuadrado: 2021/04/04 |
  |_____________________________________|

This script implements an ECG preprocessing routine (FilteringECG function),
performing the following steps:

%% Low-pass filtering (and high-pass filtering) .
Bandpass filtering is applied from 0.5 up to 40 Hz to preserve the ECG
spectral band (cascading two 5–th order Butterworth highpass and lowpass
filters).

%% Powerline interference filtering.
Powerline interference reduction is addressed by linear notch filtering
(50 Hz).

%% Baseline wander and offset removal.
Baseline wander removal is accomplished using cubic splines interpolation:
A non–overlapped sliding window of 1.2–s length is used to determine the 
knots (median of the window) from which the drift component is estimated;
baselinewander is cancelled by subtracting the estimated drift).

%% Recording spikes (artifacts of great amplitude) removal.
I am not sure how to eliminate recording peaks. I had considered using a 
threshold-based method or a median filter, but I think neither is optimal
for any ECG signal. 


In addition, this code detects if any of the channels are disconnected. 
For this, a method based on calculating the standard deviation of the 
signal derivative (or also simply the signal) is used.


%}
clear all
close all

%% LOAD AND VISUALIZATION

signal=load('ecgConditioningExample.mat').ecg;
fs=load('ecgConditioningExample.mat').fs;

raw_figure=figure(); clf;
for i=1:size(signal,2)
    subplot(size(signal,2),1,i)
    time_vector=0:1/fs:size(signal(:,i),1)/fs-1/fs;
    plot(time_vector, signal(:,i)), hold on, axis tight, grid on,
    xlabel('Time (s)'); ylabel('Amplitude'); title(['Raw ECG Channel ' num2str(i) ],'FontSize',14)
end
raw = findobj(raw_figure, 'type', 'axes', 'tag', '' );
linkaxes(raw,'x');


%% ECG SIGNAL PRE-PROCESSING 

figure_filtered=figure();clf;
for i=1:size(signal,2)
    subplot(size(signal,2),1,i)
    time_vector=0:1/fs:size(signal(:,i),1)/fs-1/fs;
    filtered_signal(:,i)= FilteringECG(signal(:,i),fs);
    plot(time_vector, filtered_signal(:,i)), hold on, axis tight, grid on,
    xlabel('Time (s)'); ylabel('Amplitude'); title(['Filtered ECG Channel ' num2str(i) ],'FontSize',14)
end
all_ha = findobj(figure_filtered, 'type', 'axes', 'tag', '' );
linkaxes(all_ha,'x');


%% DETECT DISCONNECTED CHANNELS 
for j=1:size(filtered_signal,2)
 if std(diff(filtered_signal(:,j)))< 0.01 %std(filtered_signal(:,i))< 0.01
    display(['ECG channel ' num2str(j) ' is disconnected '])
 end
end




function filtered_ecg = FilteringECG(ecg,fs)
% Time
N=length(ecg);
t=(0:N-1)/fs; 


%% Baseline wander and offset removal.
%Cubic splines interpolation, non–overlapped sliding window of 1.2–s
tw = 1.2; %window length
ecg_det=detrendSpline(ecg,fs,tw,t); 


%% Powerline interference filtering.
%Notch filtering, cut-off frequency 50 Hz.
fcancel=50; % Notch filter
ecg_notch = miFiltroNotchIIR(ecg_det,fcancel,fs,t);

%% Low-pass and High-pass filtering.
% Low-pass IIR filter cut-off frequency 40 Hz.
ecg_filtlow = miFiltroIIR(ecg_notch,fs,t);

% High-pass IIR filter cut-off frequency 0.5 Hz.
ecg_filthigh = miFiltroIIRHigh(ecg_filtlow,fs,t);

% Filter to remove possible filter trends
filtered_ecg=detrendSpline(ecg_filthigh,fs,tw,t);


end



% ### Baseline filter ###
function s_det = detrendSpline(s,fs,TBWcancell,~)

L_s=length(s);
s = s(:);
T_w = TBWcancell;
L_w=round(fs*T_w); 
t1 = 1/fs*(0:L_s-1);
t1 = t1';
t_m=[];
s_m=[];

%Median filter
for n_t=1:round(L_w/2):L_s-L_w
    t_m=[t_m t1(n_t+L_w/2-1)];
    s_m=[s_m median(s(n_t:n_t+L_w-1))];
end

%Interpolation spline
pp=csaps(t_m,s_m);
s_pp=ppval(pp,t1);
s_det=s-s_pp;
residue = s-s_det;
end

% ### IIR Notch filter ###

function ecg_notchiir = miFiltroNotchIIR(ecg,fcancel,fs,t)
Q=80;
wo = fcancel/(fs/2); %Notch frecuency
bw = wo/Q;
[b,a]=iirnotch(wo,bw); % Notch filter implementation

[H,f] = freqz(b,a,[],fs);

ecg_notchiir=filter(b,a,ecg);
residue = ecg_notchiir-ecg; 
end

% ### Low-pass IIR filter ###
function ecg_filtiir = miFiltroIIR(ecg,fs,~)
fc = 40; 
orden = 5;

[b,a] = butter(orden,fc/(fs/2));
[~,~] = freqz(b,a,[],fs);

ecg_filtiir = filtfilt(b,a,ecg);
residue = ecg_filtiir-ecg;
end


% ### High-pass IIR filter ###
function ecg_filthigh = miFiltroIIRHigh(ecg,fs,~)
fc = 0.5; 
orden = 5;

[b,a] = butter(orden,fc/(fs/2),'high');
[~,~] = freqz(b,a,[],fs);

ecg_filthigh = filtfilt(b,a,ecg);
residue = ecg_filthigh-ecg;
end

% ### High-amplitude noise filter ###
%function high_amp_noise = Recording_spikes(s,fs,thres)


%end

