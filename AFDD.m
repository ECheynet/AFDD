function [phi,fn,zeta,varargout] = AFDD(Az,t,Nmodes,varargin)
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,varargin) calculate the mode shapes, 
% eigen frequencies and modal damping ratio of the acceleration data using
% the Automated Frequency Domain Decomposition (AFDD) method
% which is based on the Frequency Domain Decomposition (FDD) [1,2]
% 
%% Input
%  * Az: acceleration data. Matrix of size [Nyy x N] where Nyy
% is the number of sensors, and N is the number of time steps
%  * fs: sampling frequencies
%  * fn: Vecteur "target eigen frequencies". ex: fn = [f1,f2,f3]
% Optional inputs as inputParser:
%  * M: [1 x 1 ]  integer.  number of FFT points
%  * PickingMethod: automated or manual peak picking ('auto' or 'manual')
%  * fn [1 x Nmodes]:  prescribed or not prescribed eigen frequencies (empty or scalar/vector)
%  * ModeNormalization: [1 x 1 ]: 1 or 0: option for mode normalization 
%  * dataPlot: 1 or 0: option for intermediate data plot (e.g. checking procedure)
%  * Ts: [1 x 1]: float: option for duration of autocorrelation function (for estimation of damping ratio only)
%  * LL: [1 x Nmodes]: float: option for selectin the lower boundary for cut-off frequency
%  * UL: [1 x Nmodes]: float: option for selectin the uper boundaty for cut-off frequency
% 
%% Output
% phi: matrix of the measured mode shapes. Size(phi)=  [Nyy x numel(fn)]
% fn: matrix of the measured eigen frequencies. Size(phi)=  [Nyy x numel(fn)]
% phi: matrix of the measured mode shapes. Size(phi)=  [Nyy x numel(fn)]
% 
%% Example
% See Example2.m joined in the FileExchange submission
% 
%% Syntax
% [phi,fn,zeta] = AFDD(Az,t,Nmodes)
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,'M',512) -> can be 512 or any integer
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,'PickingMethod','manual')
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,'fn',[0.5]) ->[0.5] is an example
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,'ModeNormalization',1)
% [phi,fn,zeta] = AFDD(Az,t,Nmodes,'dataPlot',0)
% 
%% Author Info
% E. Cheynet - UiB - last modified 17-11-2020
% 
%% References
% [1] Brincker, R.; Zhang, L.; Andersen, P. (2001). 
% "Modal identification of output-only systems using frequency domain 
% decomposition". Smart Materials and Structures 10 (3): 441. 
% doi:10.1088/0964-1726/10/3/303.
% 
% [2] Brincker, R., Zhang, L., & Andersen, P. (2000, February).
% Modal identification from ambient responses using frequency domain
% decomposition. In Proc. of the 18*選nternational Modal Analysis Conference
% (IMAC), San Antonio, Texas.

%% Input parser for options

p = inputParser();
p.CaseSensitive = false;
p.addOptional('M',2^(nextpow2(numel(t)/8))) % number of FFT points
p.addOptional('PickingMethod','auto'); % automated or manual peak picking ('auto or manual')
p.addOptional('UB',[]); % upper boundary for cut-off frequency
p.addOptional('LB',[]); % lower boundary for cut-off frequency
p.addOptional('fn',[]) ; % prescribed or not prescribed eigen frequencies (empty or scalar/vector)
p.addOptional('ModeNormalization',1) % option for mode normalization (1 or 0)
p.addOptional('dataPlot',0) % option for intermediate data plot (e.g. checking procedure) (1 or 0)
p.addOptional('Ts',30) % option for duration of autocorrelation function (for estimation of damping ratio only)
p.parse(varargin{:});

% shorthen the variables name
M = p.Results.M ;
fnTarget = p.Results.fn ;
PickingMethod = p.Results.PickingMethod;
ModeNormalization = p.Results.ModeNormalization;
dataPlot = p.Results.dataPlot;
Ts = p.Results.Ts;
LB = p.Results.LB;
UB = p.Results.UB;
%% Check error and unexpected inputs
% check if M is empty
if isempty(M),
    warning('M is specified as empty. The default value is taken instead');
    M = 2^(nextpow2(numel(t)/8));
end
% Check picking method:
if ~strcmpi(PickingMethod,'auto') && ~strcmpi(PickingMethod,'manual')
    error(' ''PickingMethod'' is not recognized. It must be a either ''auto'' or ''manual'' ')
end
% Check if numel(fn) is different from Nmodes
if ~isempty(fnTarget) && numel(fnTarget) ~=Nmodes,
    error(' The number of eigen frequencies specified fn is different from Nmodes. They must be identical');
end
% Check if dataPLot is different from 0 or 1
if dataPlot~=1 && dataPlot~=0
    error(' The value of ''dataPlot'' is not recognized. it must be 0 or 1 ');
end

if isempty(LB),
    LB = nan(1,Nmodes);
elseif numel(LB)<Nmodes,
    error('numel(LB) ~= Nmodes')
elseif min(size(UB))>1,
    error('LB should be a vector, not a matrix')
end
if isempty(UB),
    UB = nan(1,Nmodes);
elseif numel(UB)<Nmodes,
    error('numel(UB) ~= Nmodes')
elseif min(size(UB))>1,
    error('UB should be a vector, not a matrix')
end

if ~issorted(fnTarget),    warning('fnTarget are now sorted'); end

%% Pre-processing
[Nyy,N] = size(Az);
fs = 1/median(diff(t));
z = linspace(0,1,Nyy); % normalized span
%% Computation of the spectral matrix G
%  size(G) is [N x Nyy x Nyy]
if rem(M,2),
    G = zeros(Nyy,Nyy,round(M/2));
else
    G = zeros(Nyy,Nyy,round(M/2)+1);
end
for ii=1:Nyy,
    for jj=1:Nyy,
        [G(ii,jj,:),f] = cpsd(Az(ii,:),Az(jj,:),M,round(M/2),M,fs);
    end
end
%% Application of SVD to G
U =zeros(size(G));
S =zeros(Nyy,size(G,3));
V =zeros(size(G));
for ii=1:size(G,3),
    [U(:,:,ii),diagMat,V(:,:,ii)] = svd(G(:,:,ii));
    S(:,ii) = diag(diagMat);
end

if dataPlot==1,
    figure
    loglog(f,S)
    ylim([1e-10,max(S(1,:))*10]);
    xlabel('Frequency (Hz)')
    ylabel('Singular values of the PSD matrix')
end


%% interpolation to improve accuracy of peak picking and damping estimation
Ninterp=3;
newF = linspace(f(1),f(end),Ninterp*numel(f));
newS = interp1(f,S(1,:),newF,'pchip');
newS = newS./max(newS); % normalized power spectral density
newU = interp2(z,f(:),squeeze(U(:,1,:))',z,newF(:))';

if nargout==4,
    Sfdd.S = newS;
    Sfdd.f = newF;
    varargout{1} = Sfdd;
end


%% Peak picking algorithm
if isempty(fnTarget),
    if strcmpi(PickingMethod,'auto'),
        indMax = pickpeaks(newS,Nmodes,0);
        if dataPlot==1,
            figure
            loglog(newF,newS,'k-',newF(indMax),newS(indMax),'ko','markerfacecolor','r')
            ylim([min(newS(1,:)),max(newS(1,:))*10]);
            xlabel('Frequency (Hz)')
            ylabel('Normalized PSD ')
            legend('1st Singular values','peaks','location','NorthWest')
        end
    elseif strcmpi(PickingMethod,'manual'),
        indMax = manualPickPeaking(newF,newS,Nmodes);
    else
        error('pick-peaking method unknown. Please choose between ''auto'' and ''manual'' ');
    end
    % Eigen-frequencies
    fn = newF(indMax);
    % Mode shapes
    phi = zeros(numel(indMax),Nyy);
    for ii=1:numel(indMax),
        phi(ii,:) = real(newU(:,indMax(ii)));
    end
else
    fn = fnTarget;
    % Mode shapes
    phi = zeros(Nmodes,Nyy);
    for ii=1:Nmodes,
        [~,indF]=min(abs(newF-fn(ii)));
        phi(ii,:) = real(newU(:,indF));
    end
end
% normalisation of the modes
if ModeNormalization==1,
    for ii=1:size(phi,1),
        phi(ii,:) = phi(ii,:)./max(abs(phi(ii,:)));
    end
end


%% Get damping ratio
zeta = zeros(size(fn));


% sort eigen frequencies
[fn,indSort]=sort(fn);
phi = phi(indSort,:);

for ii=1:numel(fn),
    [x] = generate_response_SDOF(newF,t,newS,fn(ii),LB(ii),UB(ii));
    % We want segments of 30 seconds for the autocorrelation function
    method = 1; % cross-covariance calculated with ifft
    [IRF,newT] = NExT(x,median(diff(t)),Ts,method);
    % get the envelop of the curve with the hilbert transform:
    envelop = abs(hilbert(IRF));    envelop(1)=IRF(1);
    wn = 2*pi*fn(ii); % -> obtained with peak picking method (fast way)
    [zeta(ii)] = expoFit( envelop,newT,wn);
end


    
    function [x] = generate_response_SDOF(freq,t,S,fn,LB,UB)
% function [x] = generate_response_SDOF(freq,t,S,fc,LB,UB) generate 
% the time damain history of a SDOF system with an eigen frequency fc and
% a spectral domain limited between LB and UB. The atrget spectrum is S,
% associated with the frequency freq.


        % lower boundary for selected peak
        if isnan(LB),
            f_lower = fn*0.9;
        else
            f_lower = LB;
        end
        
        % lower boundary for selected peak
        if isnan(UB),
            f_upper = fn*1.1;
        else
            f_upper = UB;
        end
        
        [~,indLB]=min(abs(freq-f_lower)); % find left boundary
        [~,indUB]=min(abs(freq-f_upper)); % find right boundary
        % Time series generation - Monte Carlo simulation
        Nfreq = numel(S(indLB:indUB));
        df = median(diff(freq));
        w = 2*pi.*freq(:);
        A = sqrt(2.*S(indLB:indUB).*df);
        B =cos(w(indLB:indUB)*t + 2*pi.*repmat(rand(Nfreq,1),[1,N]));
        x = A*B;
    end
    function [Fp] = manualPickPeaking(f,S,Nmodes)
        % original author:  Mohammad Farshchin
        % FileExchange submission: https://se.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-/content/FDD.m
        %%
        display('Peak selection procedure')
        display('a: Draw rectangles around peaks while holding left click')
        display('b: Press "Space" key to continue the peak selection')
        display('c: Press "any other key" if you have selected a peak by mistake and want to ignore it')
        
        clf;close all;
        figure('units','normalized','outerposition',[0,0,1,1])
        plot(f,mag2db(S))
        grid on
        ylim([min(mag2db(S)),max(mag2db(10*S))])
        xlim([f(2),f(end)])
        hold on
        xlabel('Frequency (Hz)')
        ylabel('1st Singular values of the PSD matrix (db)')
        Fp=[];% Frequencies related to selected peaks
        while numel(Fp)<Nmodes
            myRec=drawrectangle;                                                                          % Draw a rectangle around the peak
            [~,P1]=min(abs(f-myRec.Position(1)));
            [~,P2]=min(abs(f-(myRec.Position(1)+myRec.Position(3))));
            [~,P3]=max(S(P1:P2));
            indPeak=P3+P1-1;                                                                         % Frequency at the selected peak
            scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','b','MarkerFaceColor','b')         % Mark this peak
            pause;
            key=get(gcf,'CurrentKey');
            if strcmp(key,'space'),
                % Press space to continue peak selection
                Fp=[Fp,indPeak];
                scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','g','MarkerFaceColor','g')      % Mark this peak as green
            else
                % Press any other key to ignore this peak
                scatter(f(indPeak),mag2db(S(indPeak)),'MarkerEdgeColor','r','MarkerFaceColor','r')      % Mark this peak as red
            end
        end
        % Number selected peaks, respectively
        Fp=sort(Fp);
        pause(0.01);
    end
    function [IRF,t] = NExT(y,dt,Ts,method)
        % [IRF] = NExT(y,ys,T,dt) implements the Natural Excitation Technique to
        % retrieve the Impulse Response FUnction (IRF) from the cross-correlation
        % of the measured output y.
        %
        % [IRF] = NExT(y,dt,Ts,1) calculate the IRF with cross-correlation
        % calculated by using the inverse fast fourier transform of the
        % cross-spectral power densities  (method = 1).
        %
        % [IRF] = NExT(y,dt,Ts,2) calculate the IRF with cross-correlation
        % calculated by using the unbiased cross-covariance function (method = 2)
        %
        %
        % y: time series of ambient vibrations: vector of size [1xN]
        % dt : Time step
        % method: 1 or 2 for the computation of cross-correlation functions
        % T: Duration of subsegments (T<dt*(numel(y)-1))
        % IRF: impusle response function
        % t: time vector asociated to IRF
        %%
        if nargin<4, method = 2; end % the fastest method is the default method
        if ~ismatrix(y), error('Error: y must be a vector or a matrix'),end
        [Ny,N1]=size(y);
        if Ny>N1,
            y=y';
            [Ny,N1]=size(y);
        end
        % get the maximal segment length fixed by T
        M1 = round(Ts/dt);
        switch method
            case 1
                clear IRF
                for pp=1:Ny,
                    for qq=1:Ny,
                        y1 = fft(y(pp,:));
                        y2 = fft(y(qq,:));
                        h0 = ifft(y1.*conj(y2));
                        IRF(pp,qq,:) = h0(1:M1);
                    end
                end
                % get time vector t associated to the IRF
                t = linspace(0,dt.*(size(IRF,3)-1),size(IRF,3));
                if Ny==1,
                    IRF = squeeze(IRF)'; % if Ny=1
                end
            case 2
                IRF = zeros(Ny,Ny,M1+1);
                for pp=1:Ny,
                    for qq=1:Ny,
                        [dummy,lag]=xcov(y(pp,:),y(qq,:),M1,'unbiased');
                        IRF(pp,qq,:) = dummy(end-round(numel(dummy)/2)+1:end);
                    end
                end
                if Ny==1,
                    IRF = squeeze(IRF)'; % if Ny=1
                end
                % get time vector t associated to the IRF
                t = dt.*lag(end-round(numel(lag)/2)+1:end);
        end
        % normalize the IRF
        if Ny==1,
            IRF = IRF./IRF(1);
        else
        end
    end
    function [zeta] = expoFit(y,t,wn)
        % [zeta] = expoFit(y,t,wn) returns the damping ratio calcualted by fiting
        % an exponential decay to the envelop of the Impulse Response Function.
        %
        % y: envelop of the IRF: vector of size [1 x N]
        % t: time vector [ 1 x N]
        % wn: target eigen frequencies (rad/Hz) :  [1 x 1]
        % zeta: modal damping ratio:  [1 x 1]
        %  optionPlot: 1 to plot the fitted function, and 0 not to plot it.
        %%
        
        % Initialisation
        guess = [1,1e-2];
        % simple exponentiald ecay function
        myFun = @(a,x) a(1).*exp(-a(2).*x);
        % application of nlinfit function
        assert(license('test','Statistics_Toolbox')==1,'The function expoFit requires Matlab Statistics Toolbox.')
        coeff = nlinfit(t,y,myFun,guess);
        % modal damping ratio:
        zeta = abs(coeff(2)./wn);
    end
end