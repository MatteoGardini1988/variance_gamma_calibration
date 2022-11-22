clc
clear variables
close all

pathImages = '../../images';
extension = 'epsc';

% For reproducibility
SEED = 4;
rng(SEED);

% Load Forward Germany prices
load('GermanyFwdPrices.mat')

% Get prices
Date = datetime(2019,11,19);
idx = germanyFwdPrices.Data <= Date;
dataGermayForward = germanyFwdPrices.DEBY2021(idx);

data = diff(log(dataGermayForward));

% Remove where the increment is zero (this is importan to avoid errors)
data(data==0) = [];

dt = 1/252;
M = mean(data);
V = var(data);
S = skewness(data);
K = kurtosis(data);
sigma = sqrt(V/dt);
nu = (K/3 -1)*dt;
theta = (S* sigma * sqrt(dt))/(3* nu );

%% VG MLE
pdf_VG = @(data,theta,nu,sigma)VGdensity_2(data,theta,nu,sigma,dt);
%%

x = linspace(min(data),max(data),500);
fx_st = VGdensity_2(x,theta,nu,sigma,dt);

% haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
% histogram(data,'Normalization','pdf');
% hold on
% plot(x,fx_st)

%% Optimize the Likelihood
% start = [theta,nu,sigma,mu];
% lb = [-intmax 0 0 -intmax ];
% ub = [intmax intmax intmax intmax ];
lb = [-intmax 0 0];
ub = [intmax intmax intmax];
options = statset ('MaxIter' ,50000,'MaxFunEvals',50000);

% Set the starting point
start = [theta,nu,sigma];

% Calibation
params = mle(data,'pdf',pdf_VG , ...
    'start',start,'lower',lb ,...
    'upper',ub,'options',options );

theta_hat = params(1);
nu_hat = params(2);
sigma_hat = params(3);


save('hist_fitting.mat','theta_hat',...
    'nu_hat','sigma_hat');

%% Compute the density
fx = VGdensity_2(x,theta_hat,nu_hat,sigma_hat,dt);

haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
h2 = histogram(data,'Normalization','pdf');
hold on
plot(x,fx,'-','Color',[0 0.5 0],'LineWidth',3)
plot(x,fx_st,'-','Color',[0.8 0 00],'LineWidth',3)
h2.FaceColor = [0 0.5 0];
h2.FaceAlpha = 0.1;
legend('Market data','Fitted \Theta','Starting \Theta_{0}');
set(gca,'FontSize',20);
xlim([-0.06 0.06]);
%saveas(haic,fullfile(pathImages,'FittingHistorical'),extension);

%% Plot
nSim = 500;
nDates = numel(dataGermayForward);
dt = 1/252;

params_sim(1) = theta_hat;
params_sim(2) = nu_hat;
params_sim(3) = sigma_hat;

X = VG_simulation(nSim,nDates,dt,params_sim);

xLOC = X(:,2:end);
Y = diff(xLOC(:));

figure
qqplot(data,Y);

S = dataGermayForward(1).*exp(X);

figure
plot(S');
hold on
plot(dataGermayForward,'LineWidth',3,'Color',[0 0 0]);

