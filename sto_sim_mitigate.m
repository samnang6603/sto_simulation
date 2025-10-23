% STO_estimation.m
clear,
rng('shuffle')
nSTOs = [-3 -3 2 2]; % Number of samples corresponding to STO
CFOs = [0 0.5 0 0.5]; SNRdB=40; MaxIter=10; % CFOs, SNR, # of iteration
Nfft=128; Ng=Nfft/4; % FFT size and GI (CP) length
Nofdm=Nfft+Ng; % OFDM symbol length
Nbps=2; M=2^Nbps; % Number of bits per (modulated) symbol
%mod_object = modem.qammod(’M’,M, ’SymbolOrder’,’gray’);
Es=1; A=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
N=Nfft; com_delay=Nofdm/2; Nsym=100;
%rand(’seed’,1); randn(’seed’,1);

for i=1:length(nSTOs)
    nSTO=nSTOs(i); CFO=CFOs(i);
    x = []; % Initialize a block of OFDM signals
    for m=1:Nsym % Transmit OFDM signals
        %msgint=randint(1,N,M);
        msgint = randi([0 M-1],N,1);
        %Xf = A*modulate(mod_object,msgint);
        Xf = A*qammod(msgint,M).';
        xt = ifft(Xf,Nfft); 
        x_sym = [xt(end-Ng+1:end) xt]; % IFFT & Add CP
        x = [x x_sym];
    end
    y = x; % No channel effect
    y_CFO = add_CFO(y,CFO,Nfft); 
    y_CFO = y;
    y_CFO_STO = add_STO(y_CFO,-nSTO);
    Mag_cor= 0; 
    Mag_dif= 0;
    for iter=1:MaxIter
        y_aw = awgn(y_CFO_STO,SNRdB,'measured'); % AWGN added
        % Symbol Timing Acqusition
        [STO_cor,mag_cor]=STO_by_correlation(y_aw,Nfft,Ng,com_delay);
        [STO_dif,mag_dif]=STO_by_difference(y_aw,Nfft,Ng,com_delay);
        Mag_cor= Mag_cor+mag_cor; Mag_dif= Mag_dif+mag_dif;
    end % End of for loop of iter
    [Mag_cor_max,ind_max] = max(Mag_cor); nc= ind_max-1-com_delay;
    [Mag_dif_min,ind_min] = min(Mag_dif); nd= ind_min-1-com_delay;
    nn=-Nofdm/2+(0:length(Mag_cor)-1);
    subplot(220+i),
    plot(nn,Mag_cor,'b:', nn,Mag_dif,'r'), hold on
    str1=sprintf('Cor(b-)/Dif(r:) for nSTO=%d, CFO=%1.2f',nSTO,CFO);
    title(str1);
    xlabel('Sample'),
    ylabel('Magnitude');
    stem(nc,Mag_cor(nc+com_delay+1),'b') % Estimated STO from correlation
    stem(nd,Mag_dif(nd+com_delay+1),'r') % Estimated STO from difference
    stem(nSTO,Mag_dif(nSTO+com_delay+1),'k.') % True STO
end % End of for loop of i

function y_CFO=add_CFO(y,CFO,Nfft)
% add CFO (carrier frequency offset)
% y : Received signal
% CFO = IFO (integral CFO) + FFO (fractional CFO)
% Nfft = FFT size
nn=0:length(y)-1; y_CFO = y.*exp(1j*2*pi*CFO*nn/Nfft); % Eq.(5.7)
end


function y_STO=add_STO(y, nSTO)
% add STO (symbol time offset)
% y : Received signal
% nSTO : Number of samples corresponding to STO
if nSTO>=0
    y_STO=[y(nSTO+1:end) zeros(1,nSTO)]; % advance
else
    y_STO=[zeros(1,-nSTO) y(1:end+nSTO)]; % delay
end
end

function [STO_est, Mag]=STO_by_correlation(y,Nfft,Ng,com_delay)
% estimates STO by maximizing the correlation between CP (cyclic prefix)
% and rear part of OFDM symbol
% Input: y = Received OFDM signal including CP
% Ng = Number of samples in Guard Interval (CP)
% com_delay = Common delay
% Output: STO_est = STO estimate
% Mag = Correlation function trajectory varying with time
Nofdm=Nfft+Ng; % OFDM symbol length
if nargin<4
    com_delay = Nofdm/2; 
end
nn=0:Ng-1;
% Take correlation
%{

  Ng-1+del
    --
    \  |y[n+i]*conj(y[n + N + i])| = T
    /_
  i = del

  Max(T)

del = com_delay
Ng-1+del = length(yy)
N = Nfft
%}
yy = y(nn+com_delay)*y(nn+com_delay+Nfft)'; % Correlation
maximum=abs(yy);

for n=1:Nofdm
    n1 = n-1;
    yy1 = y(n1+com_delay)*y(n1+com_delay+Nfft)';
    yy2 = y(n1+com_delay+Ng)*y(n1+com_delay+Nfft+Ng)';
    yy = yy-yy1+yy2; Mag(n)=abs(yy); % Eq.(5.12)
    if Mag(n)>maximum
        maximum=Mag(n);
        STO_est=Nofdm-com_delay-n1;
    end
end
end

function [STO_est,Mag]=STO_by_difference(y,Nfft,Ng,com_delay)
% estimates STO by minimizing the difference between CP (cyclic prefix)
% and rear part of OFDM symbol
% Input: y = Received OFDM signal including CP
% Ng = Number of samples in CP (Guard Interval)
% com_delay = Common delay
% Output: STO_est = STO estimate
% Mag = Correlation function trajectory varying with time
Nofdm=Nfft+Ng; minimum=100; STO_est=0;
if nargin<4
    com_delay = Nofdm/2; 
end
for n=1:Nofdm
    nn = n+com_delay+(0:Ng-1); tmp0 = abs(y(nn))-abs(y(nn+Nfft));
    Mag(n) = tmp0*tmp0'; % Squared difference by Eq.(5.11)
    if Mag(n) 
        minimum=Mag(n); 
        STO_est=Nofdm-com_delay-(n-1); 
    end
end
end

