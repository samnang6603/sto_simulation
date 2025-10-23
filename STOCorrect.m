function [y,timingError] = STOCorrect(x,sps,zeta,nbw,Kp,TED)

if nargin < 3
    [zeta,nbw,Kp,TED] = deal(1,0.01,2.7,'Zero-Crossing');
elseif nargin < 4
    [nbw,Kp,TED] =  deal(0.01,2.7,'Zero-Crossing');
elseif nargin < 5
    [Kp,TED] = deal(2.7,'Zero-Crossing');
else
    TED = 'Zero-Crossing';
end

% Constant to tune interpolation filter
alpha = 0.5;
expFactor = [11 10];

% PLL param
inputFrameLen = length(x);
maxOutputFrameLen = ceil(inputFrameLen*expFactor(1))/(expFactor(2)*sps);

% interpolation filter
hInterp = [0,0,1,0; % constant
    -alpha, 1+alpha, -(1-alpha), -alpha;  % Linear
     alpha,  -alpha,     -alpha,  alpha]; % Quadratic
mu = 0;
interpFiltState = zeros(3,1,'like',1i);

% Additional loop filter param
% C.56 and C.57 page 736 Rice's book
K0 = -1; % modulo-1
theta = nbw/(sps*(zeta + 1/(4*zeta)));
den = Kp*K0*(1 + 2*zeta*theta + theta^2);
K1 = 4*zeta*theta/den; % Proportional gain
K2 = 4*theta^2/den; % Integrator gain

% Symbol Sync
strobe = 0;
numStrobe = 0;
overflowFlag = false;
strobeHistory = false(1,sps);
bufferTED = zeros(1,sps,'like',1i);
symHolder = zeros(maxOutputFrameLen,1,'like',1i);
symHolderLen = maxOutputFrameLen;
timingError = zeros(inputFrameLen,1);
loopFiltState = 0;
ncocnt = 0;

for sampIdx = 1:inputFrameLen

    % if symHolder buffer is filled with strobes (aligned syms), then the
    % DPLL has gone rogue and wants to keep sampling beyond what it's
    % supposed to
    if (numStrobe == symHolderLen) && strobe
        overflowFlag = true;
        break;
    end

    numStrobe = numStrobe + strobe;
    timingError(sampIdx) = mu;

    % Figure 8.4.2 on page 449 in Rice's book [1].
    [intFiltOut,interpFiltState] = interpFilt(x(sampIdx),hInterp,mu,interpFiltState);

    if strobe
        symHolder(numStrobe) = intFiltOut;
    end

    if numStrobe > symHolderLen
        overflowFlag = true;
        break;
    end

    e = elTED(intFiltOut,bufferTED,strobe,strobeHistory);

    switch sum([strobeHistory(2:end), strobe])
        case 0
            % NOP
            % Skip current sample if no strobe
        case 1
            % Shift TED buffer regularly if ONE strobe across N samples, i.e.,
            bufferTED = [bufferTED(2:end), intFiltOut];
        otherwise
            % Stuff a missing sample if TWO strobes across N samples, i.e.,
            % strobeHistory(2:end) = [1, 0, ..., 0] & strobe = 1
            bufferTED = [bufferTED(3:end), 0, intFiltOut];
    end

    % loop filter
    [v,loopFiltState] = loopFilt(e,loopFiltState,K1,K2);

    % interpolation control
    [ncocnt,mu,strobe,strobeHistory] = interpControl(v,sps,ncocnt,...
                                        mu,strobe,strobeHistory);

end
y = symHolder(1:numStrobe,1);
if overflowFlag
    warning('Overflow: Dropping Synchronization')
end




% %% TED
% switch TED
%     case 'Zero-Crossing'
% 
%     case 'Gardner'
% 
%     case 'Early-Late'
% 
%     case 'Mueller-Muller'
% 
%     otherwise
% end

end

%% Local Fcn
function timeerr = zcTED()



end

function timeerr = gdnTED()

end

function timeerr = elTED(x,TEDbuff,strobe,strobeHistory)
% TED calculation occurs one sample after a strobe
    if (~strobe) && strobeHistory(end) && all(~strobeHistory(1:end-1))
      timeerr = real(TEDbuff(end)) * (real(x) - real(TEDbuff(end-1))) + ...
        imag(TEDbuff(end)) * (imag(x) - imag(TEDbuff(end-1)));
    else
      timeerr = 0;
    end
end

function timeerr = mmTED()

end


function [y,filtState] = interpFilt(x,h,mu,filtState)
% Piecewise parabolic interpolator in farrow structure with alpha =
% 0.5. Refer to (8.72)-(8.73) on page 468 and Figure 8.4.16 on page 471
% in Rice's book [1].
xSeq = [x; filtState];
a = h*xSeq;
b = [1; mu; mu^2];
y = sum(a.*b);
filtState = xSeq(1:3);
end

function [ncocnt,mu,strobe,strobeHistory] = interpControl(v,sps,ncocnt,...
                                            mu,strobe,strobeHistory)
% Modulo-1 counter interpolation controller which generates/updates
% strobe signal (strobe) and fractional interpolation interval
% (mu). Refer to Section 8.4.3 and Figure 8.4.19 in Rice's book [1].

W = v + 1/sps;
strobeHistory = [strobeHistory(2:end), strobe]; % add in new strobe sample
strobe = ncocnt < W;
if strobe
    mu = ncocnt/W;
end
ncocnt = mod(ncocnt - W,1); % update counter

end

function [v,filtState] = loopFilt(timeerr,filtState,K1,K2)
integFiltOut = timeerr*K2 + filtState;
v = timeerr*K1 + integFiltOut;
filtState = integFiltOut;
end