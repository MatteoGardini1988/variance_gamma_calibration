clear variables
clc
close all

% RISK-NEUTRAL CALIBRATION and 

% extension = 'eps';
extension = 'eps';

seed = 4;
rng(seed);

pathImages = '../images';

% Load data
load('PwC_OptData.mat');
load('PwC_FwdData.mat');

% Make Categorical
dati_opt.Underlying = categorical(dati_opt.Underlying);
dati_opt.Type = categorical(dati_opt.Type);
dati_fwd.Contract = categorical(dati_fwd.Contract);



% Take the calendar Fwd
idxSel = dati_fwd.Contract == 'DEBY';
dati_fwd = dati_fwd(idxSel,:);

% Take the 2021
idxSel = strcmp(dati_fwd.DeliveryPeriod,'2021.01');
dati_fwd = dati_fwd(idxSel,:);

% Price of Forward
F0 = dati_fwd.SettlementPrice;

% Take option on Calendadar
idxSel = dati_opt.Underlying == 'DEBY 2021.01';
dati_opt = dati_opt(idxSel,:);

% Take only Call options
idxSel = dati_opt.Type == 'C';
dati_opt = dati_opt(idxSel,:);

% Take only ITM Call options
% idxSel = dati_opt.Strike >= F0;
% dati_opt = dati_opt(idxSel,:);

r = 0.01;
t0 = datetime(dati_fwd.TradingDate,'InputFormat','yyyy-MM-dd');
maturities = datetime(dati_opt.ExpiryDate,'InputFormat','yyyy-MM-dd'); 
T  = yearfrac(t0,maturities,13);
K = dati_opt.Strike;
P = dati_opt.SettlementPrice;

Tunique  = unique(T);

idxSel = T == T(end);

haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
sigmaImp = blkimpv(F0,K(idxSel),r,T(idxSel),P(idxSel));
toPlot = [K(idxSel) sigmaImp];
toPlot = sortrows(toPlot,1);
plot(toPlot(:,1),toPlot(:,2),'-x','LineWidth',1.5,'Color',[0 0 0]);
title('');
xlabel('Strike Price [EUR/MWh]');
ylabel('\sigma - Implied Volatility');
title('Implied Volatility Smile');
set(gca,'FontSize',20);

% saveas(haic,fullfile(pathImages,'VolatilitySmile'),extension);


%% Calibration of the Model
% Starting point
x0 = [1.05 0.02 0.2];

% Upper and lower bounds
lb = [-2.7 0.01 0.02];
ub = [2.7 0.8 1.1];


% calibration Routine
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',5000,'MaxIterations',2000);

% Perform Calibration
[params,resnorm,residual] = ...
    lsqnonlin(@(x)CalibFunction(x,F0,T,r,K,P),x0,lb,ub,options);

theta = params(1);
nu = params(2);
sigma = params(3);

save('riskneutral_fitting.mat','theta',...
    'nu','sigma');

disp('Parameters');
fprintf('theta: %.4f\n',theta);
fprintf('nu: %.4f\n',nu);
fprintf('sigma: %.4f\n',sigma);

%% Compute implied volatilities

implVol = blkimpv(F0,K,r,T,P);
dati_opt.ImplVol = implVol;

%% Market repricing

dati_opt.Repricing = nan(size(P));


omega = (1/nu).*...
    log(1-theta*nu - sigma*sigma*nu/2);

% Get distinct Maturities
Tdistinct = unique(T);
Nmaturities = length(Tdistinct);


haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
left_color = [0 0 0];
right_color = [0 0 0];
set(haic,'defaultAxesColorOrder',[left_color; right_color]);

for i = 1:Nmaturities
    % Get maturities
    idx = T == Tdistinct(i);
    
    
    % perform pricing via FFT
    [CallPrices,LogStrikes] = FFTPricing(Tdistinct(i),r,...
        @(w)phi_vg(w,F0,r,omega,Tdistinct(i),theta,nu,sigma));
    
    % Get call Prices corrisponding to given prices
    K_sel = K(idx);
    
    % Prezzo delle Call
    Call_Prices = interp1(LogStrikes,CallPrices,log(K_sel));
    
    dati_opt.Repricing(idx) = Call_Prices;
    
    errore = Call_Prices - dati_opt.SettlementPrice(idx);
    
    subplot(sqrt(Nmaturities),sqrt(Nmaturities),i);
    plot(dati_opt.Strike(idx),dati_opt.SettlementPrice(idx),'o','LineWidth',1.5,'Color',[0 0.5 0]);
    hold on
    plot(dati_opt.Strike(idx),dati_opt.Repricing(idx),'x','LineWidth',1.5,'Color',[0.8 0 0]);
    yyaxis right
    plot(dati_opt.Strike(idx),errore,'.k');
    yyaxis left
    legend('Mkt','Model','err');
    xlabel('Strike');
    ylabel('Price');
    set(gca,'FontSize',15);
    
    
    
end

% saveas(haic,fullfile(pathImages,'MktRepricing'),extension);


%% Monte Carlo pricing and volatility surface plot
Tmax = max(Tdistinct);
nSim = 1e4;
dt = 1/252;
nDates = ceil(Tmax/dt)+1;

Y = VG_simulation(nSim,nDates,dt,params);

tt = 0:dt:Tmax;

F = F0.*exp((r+omega).*repmat(tt,nSim,1) + Y);

%% Pricing
disp('Pricing...');
Ks = min(K):0.5:max(K);

nPoints = 20;
Ts = linspace(0,Tmax,nPoints);


[STRIKE,MATURITY] = meshgrid(Ks,Ts);

nMaturities = numel(Ts);
nStrikes = numel(Ks);

[CALL,IVSURFACE] = deal(nan(nMaturities,nStrikes));

for t = 1:nMaturities
    
    [~,idxT] = min(abs(tt-Ts(t)));
    
    for k = 1:nStrikes
        
        payoff = max(F(:,idxT)-Ks(k),0);
        
        CALL(t,k) = exp(-r.*Ts(t)).*mean(payoff);
        IVSURFACE(t,k) = blkimpv(F0,Ks(k),r,Ts(t),CALL(t,k));
        
    end
    
end
disp('OK Pricing...');

%% Some plots
haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
h = surf(MATURITY,STRIKE,CALL);
h.FaceAlpha = 0.5;
hold on
scatter3(T,K,P,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0.8 0 0]);
xlabel('Maturity');
ylabel('Strike price');
title('Market Repricing');
view([45 56 45])
set(gca,'FontSize',20);
% saveas(haic,fullfile(pathImages,'RepricingPlot'),extension);



haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
h= surf(MATURITY,STRIKE,IVSURFACE);
hold on
h.FaceAlpha = 0.6;
scatter3(T,K,implVol,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 0 0]);
xlabel('Maturity');
ylabel('Strike price');
zlabel('Volatility');
title('Implied Volatility Surface');
view([-40 50 45])
set(gca,'FontSize',20);
% saveas(haic,fullfile(pathImages,'surfacePlot'),extension);

%% Plot dei Log-rendimenti dello storico

load('PwCDataForwardHist.mat');
dataForwardHist.UnderlyingName = ...
    categorical(dataForwardHist.UnderlyingName);
idxSel = (dataForwardHist.UnderlyingName == 'DEBY 2020.01');
x = diff(log(dataForwardHist.SettlementPrice(idxSel)));

idxSel = (dataForwardHist.UnderlyingName == 'DEBY 2021.01');
x =[x; diff(log(dataForwardHist.SettlementPrice(idxSel)))];


muS = mean(x);
sigmaS = std(x);

xN = muS + sigmaS.*randn(numel(x),1);


% haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
% histogram(x,'Normalization','pdf')
% hold on
% histogram(xN,'Normalization','pdf')

%% Plot dei Lon-rendimenti delle simulazioni
% dF = diff(log(F),[],2);
% muN = mean(dF(1,:));
% sigmaN = std(dF(1,:)); 
% dFN = muN + sigmaN.*randn(numel(dF(1,:)),1);
% 
% haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
% histogram(dF(1,:),'Normalization','pdf')
% hold on
% histogram(dFN(:),'Normalization','pdf')



