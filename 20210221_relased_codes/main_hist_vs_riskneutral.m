clear variables
clc
close all

% extension = 'eps';
extension = 'epsc';

seed = 4;
rng(seed);

pathImages = '../../images';

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


%% Pricing using Historical Prices vs Risk-Neutral one
load('hist_fitting.mat');
load('riskneutral_fitting.mat');



%% 
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

%saveas(haic,fullfile(pathImages,'VolatilitySmile'),extension);

%% Compute implied volatilities

implVol = blkimpv(F0,K,r,T,P);
dati_opt.ImplVol = implVol;

%% Market Repricing



dati_opt.Repricing = nan(size(P));
dati_opt.RepricingH = nan(size(P));

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
    
    [CallPricesH,LogStrikes] = FFTPricing(Tdistinct(i),r,...
        @(w)phi_vg(w,F0,r,omega,Tdistinct(i),theta_hat,nu_hat,sigma_hat));
    
    % Get call Prices corrisponding to given prices
    K_sel = K(idx);
    
    % Prezzo delle Call
    Call_Prices = interp1(LogStrikes,CallPrices,log(K_sel));
    Call_PricesH = interp1(LogStrikes,CallPricesH,log(K_sel));
    
    dati_opt.Repricing(idx) = Call_Prices;
    dati_opt.RepricingH(idx) = Call_PricesH;
    
    errore = Call_Prices - dati_opt.SettlementPrice(idx);
    erroreH = Call_PricesH - dati_opt.SettlementPrice(idx);
    
    subplot(sqrt(Nmaturities),sqrt(Nmaturities),i);
    plot(dati_opt.Strike(idx),dati_opt.SettlementPrice(idx),'o','LineWidth',1.5,'Color',[0.8 0 0]);
    hold on
    plot(dati_opt.Strike(idx),dati_opt.Repricing(idx),'x','LineWidth',1.5,'Color',[0 0.5 0]);
    plot(dati_opt.Strike(idx),dati_opt.RepricingH(idx),'x','LineWidth',1.5,'Color',[0 0 0.5]);
    yyaxis right
    plot(dati_opt.Strike(idx),errore,'x','Color',[0 0.5 0]);
    hold on
    plot(dati_opt.Strike(idx),erroreH,'x','Color',[0 0 0.5]);
    ylabel('error');
    yyaxis left
    legend('Mkt','Model RN','Model H','err RN','err H');
    xlabel('Strike');
    ylabel('Price');
    set(gca,'FontSize',15);
    
    
    
end

%saveas(haic,fullfile(pathImages,'MktRepricingHvsRN'),extension);



%%  Monte Carlo Pricing(and volatility surface plot)
Tmax = max(Tdistinct);
nSim = 5e5;
dt = 1/252;
nDates = ceil(Tmax/dt)+1;

params(1) = theta;
params(2) = nu;
params(3) = sigma;
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
%saveas(haic,fullfile(pathImages,'RepricingPlot'),extension);



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
%saveas(haic,fullfile(pathImages,'surfacePlot'),extension);
