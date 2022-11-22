% Load prices
load('GermanyFwdPrices.mat')
pathImages = '../images';
extension = 'epsc';

haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
plot(germanyFwdPrices.Data,germanyFwdPrices.DEBY2021);
title('German Forward Prices');
set(gca,'FontSize',20);
saveas(haic,fullfile(pathImages,'HistoricalSeries'),extension);


% Compute the increments:
x = diff(log(germanyFwdPrices.DEBY2021));

haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
autocorr(x);
set(gca,'FontSize',20);
% saveas(haic,fullfile(pathImages,'Autocorrelation'),extension);

% QQplot of the log returns
haic = figure('Units','normalized','OuterPosition',[0 0 1 1]);
qqplot(x);
set(gca,'FontSize',20);
% saveas(haic,fullfile(pathImages,'QQplot'),extension);