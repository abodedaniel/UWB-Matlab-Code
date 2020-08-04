    %% This code models a TR-UWB signal generator and receiver
    %  by Abode Daniel
    
    % Defining variables
    close all;
    Ep = 1;                     %Amplitude determinant
    fs = 3e10;      %Hz   sampling frequency
    Fn = fs/2;      %Hz   minimum nyquist rate
    TD = 50e-9;     %seconds     distance between reference and information pulse
    TS = 300e-9;    %seconds     distance between reference pulses
    tp = 2e-9;      %seconds     pulse width
    t=-1e-9:1/fs:1.85e-6;       %signal time
    t_ref = t;
    t_pulse = t;
    a=tp/2.5;       %fudge factor
    
    %% Genrating a single pulse
    Y=(1-(4*pi.*(t.^2))/a^2) .* exp(-2*pi.*(t.^2)/a^2) / sqrt(Ep);
    plot(t,Y,'-r');
    axis([-1e-9 1e-9 -0.6 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of a single pulse - 2nd order derivated','fontSize',14)
    
    %% Creating the reference pulses delaying by TS
    binary_seq = [1 0 1 1 0 1 0];
    for i = 1:length(binary_seq)       
            t_ref = t_ref - TS;
            Y = Y + (1-(4*pi.*(t_ref.^2))/a^2) .* exp(-2*pi.*(t_ref.^2)/a^2) / sqrt(Ep);
    end
    
    %% Creating the information pulses delaying by TD and TS
    for i = 1:length(binary_seq)
        if i == 1
            t_pulse = t_pulse - TD;
            Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(Ep); 
        end
        if i > 1
            if binary_seq(i) == 1
                t_pulse = t_pulse - TS;
                Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(Ep);
            end
            if binary_seq(i) == 0
                t_pulse = t_pulse - TS;
                Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(3);
            end
        end
    end
    
    %% The generated signal and plotting it versus time
    binary_signal = Y;
    figure
    plot(t,binary_signal,'-r')
    axis([-40e-9 1.9e-6 -0.6 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of a binary coded pulse signal','fontSize',14)
    
    figure
    plot(t,binary_signal,'-r')
    axis([-3e-9 6e-8 -0.6 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot showing only the first data (1)','fontSize',14)
    
    figure
    plot(t,binary_signal,'-r')
    axis([2.9e-7 3.6e-7 -0.6 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot showing only the second data (0)','fontSize',14)
     

    %% Generating and plotting the frequency spectrum
    
    NFFY=2.^(ceil(log(length(binary_signal))/log(2)));
    FFTY=fft(binary_signal,NFFY);
    NumUniquePts=ceil((NFFY+1)/2); 
    FFTY=FFTY(1:NumUniquePts);
    MY=abs(FFTY);
    MY=MY*2;
    MY(1)=MY(1)/2;
    MY(length(MY))=MY(length(MY))/2;
    MY=MY/length(binary_signal);
    f=(0:NumUniquePts-1)*2*Fn/NFFY;
    figure
    plot(f,20*log10(MY));xlabel('Frequency (Hz)','fontSize',14); title('Plot of Power Spectrum','fontSize',14);ylabel('Power Density (dB/Hz)','fontSize',14);
    
    
    %% Modulating the baseband signal and plotting it
    fp =4.492e9;
    [modulated_signal,tm] = modulate(binary_signal,fp,fs,'amssb');
    figure
    plot(tm,modulated_signal,'-r')
    axis([-4e-8 1.9e-6 -1.2 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of the Modulated signal','fontSize',14)
    
    figure
    plot(tm,modulated_signal,'-r')
    axis([2.95e-7 3.55e-7 -1.2 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of a part of the Modulated signal','fontSize',14)
    
    %% Generating and plotting the frequency spectrum of the modulated signal
    
    NFFY=2.^(ceil(log(length(modulated_signal))/log(2)));
    FFTY=fft(modulated_signal,NFFY);
    NumUniquePts=ceil((NFFY+1)/2); 
    FFTY=FFTY(1:NumUniquePts);
    MY=abs(FFTY);
    MY=MY*2;
    MY(1)=MY(1)/2;
    MY(length(MY))=MY(length(MY))/2;
    MY=MY/length(modulated_signal);
    f=(0:NumUniquePts-1)*2*Fn/NFFY;
    figure
    plot(f,20*log10(MY));xlabel('Frequency (Hz)','fontSize',14); title('Plot of Power Spectrum of Modulated Signal','fontSize',14);ylabel('Power Density (dB/Hz)','fontSize',14);
   
    
    %% 2.1 Squaring the received signal and plotting it
    received_signal_sqr = modulated_signal.^2;
    figure
    plot(tm,received_signal_sqr,'-r')
    axis([-4e-8 1.9e-6 -0.2 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of the recieved signal after squaring','fontSize',14)
    
    figure
    plot(tm,received_signal_sqr,'-r')
    axis([2.95e-7 3.55e-7 -0.2 1.2])
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of a part of the recieved signal after squaring','fontSize',14)
    
    %% Generating and plotting the frequency spectrum of the squared received signal
    NFFY=2.^(ceil(log(length(received_signal_sqr))/log(2)));
    FFTY=fft(received_signal_sqr,NFFY);
    NumUniquePts=ceil((NFFY+1)/2); 
    FFTY=FFTY(1:NumUniquePts);
    MY=abs(FFTY);
    MY=MY*2;
    MY(1)=MY(1)/2;
    MY(length(MY))=MY(length(MY))/2;
    MY=MY/length(received_signal_sqr);
    f=(0:NumUniquePts-1)*2*Fn/NFFY;
    figure
    plot(f,20*log10(MY));xlabel('Frequency (Hz)','fontSize',14); title('Plot of Power Spectrum of Recieve signal squared','fontSize',14);ylabel('Power Density (dB/Hz)','fontSize',14);
    
   
    %% Integrating the signal over time by batches of 100 nsec and plotting
    integrated_signal = zeros(1,length(received_signal_sqr));
    start_time = min(tm);
    next_time = 0;
    while next_time < max(tm)
        next_time = next_time + 100e-9; %every 100 nsec
        index = find((start_time <= tm) & (tm < next_time));
        integrated_signal(index) = cumtrapz(received_signal_sqr(index)); %integrating every 100 nsec
        start_time = next_time;
    end
    
    figure
    plot(tm,integrated_signal,'-r')
    ylabel('Amplitude', 'fontSize',14)
    xlabel('time in seconds','fontSize',14)
    title('Plot of the integrated signal','fontSize',14)
    
    %% Designing the decision level block
    recovered_data = zeros(1,7);
    threshold = 14;
    start_time = min(tm);
    next_time = 0;
    i = 1;
    result = zeros(1,length(binary_signal));
    while next_time < max(tm)
        next_time = next_time + 100e-9; %every 100 nsec
        index = find((start_time <= tm) & (tm < next_time));
        if(max(integrated_signal(index)) > 12 && max(integrated_signal(index) > 0))
            recovered_data(i) = 1;
            result(index) = integrated_signal(index);
        else
            recovered_data(i) = 0;
        end
        next_time = next_time + 200e-9; % delaying for 200 nsec 
        start_time = next_time;
        i = i + 1;
    end
     figure
    stem(recovered_data,'-r')
    ylabel('Amplitude', 'fontSize',14)
    title('Plot of the result','fontSize',14)
    axis([0 8 -0.1 1.1])
    disp('recovered_data : ')
    disp(recovered_data)
    
    
       
    