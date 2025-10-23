M = 4;         % Modulation order for QPSK
nSym = 5000;   % Number of symbols in a packet
sps = 4;       % Samples per symbol
timingErr = 2; % Samples of timing error
snr = 15;      % Signal-to-noise ratio (dB)

txfilter = comm.RaisedCosineTransmitFilter( ...
    OutputSamplesPerSymbol=sps);
rxfilter = comm.RaisedCosineReceiveFilter( ...
    InputSamplesPerSymbol=sps, ...
    DecimationFactor=sps/2);

symbolSync = comm.SymbolSynchronizer;

data = randi([0 M-1],nSym,1);
modSig = pskmod(data,M,pi/4);

fixedDelay = dsp.Delay(timingErr);
fixedDelaySym = ceil(fixedDelay.Length/sps);

txSig = txfilter(modSig);
delaySig = fixedDelay(txSig);

rxSig = awgn(delaySig,snr,'measured');

rxSample = rxfilter(rxSig);  
scatterplot(rxSample(1001:end),2)

rxSync = symbolSync(rxSample);
scatterplot(rxSync(1001:end),2)