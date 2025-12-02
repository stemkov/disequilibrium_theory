% AUTHOR:   Michael Cortez, Florida State University
% DATE:     December 2025
% PURPOSE:  Plot time series for formulas F2-F7, assuming linearly decreasing climate

% INPUTS:   Parameters for ecosystem function formula (alpha, beta, gamma)
%           Delta0 = Initial disequilibirum when climate change begins
%           C0 = initial climate
%           r = rate of change of climate
%           lambda = response rate of community

clear all

figure('Position',[10 10 1200 600]);
textsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters (used in all panels)
tvals = 0:0.01:10;
r = 1;
lambda = 1;
alpha = 20;
gamma =1;
C0 = 5;
C = C0 - r*tvals;
Colors = {'b','c', 'k', 'm','r'};
Delta0 = -r*lambda - [-3, -2, -1, -0.5, 2];
F = zeros(length(Delta0),length(tvals));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel A: F2

%%% Generate and plot time series
subplot(2,4,1) 
hold on
beta = 2;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); set(gca,'XTick',[0 2 4 6])
ylim([0 30]); set(gca,'YTick',[0 10 20 30])
xlabel(' ')
ylabel({'Ecosystem function'}) 
title('a: F_2') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('$\Delta_0>0, \bar{Q}<0$','$\Delta_0>0, \bar{Q}<0$',...
    '$\Delta_0=0, \bar{Q}<0$','$\Delta_0<0, \bar{Q}<0$','$\Delta_0<0, \bar{Q}>0$',...
    'position',[0.35 0.1 1 1],'FontSize', 16,'Interpreter','Latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel B: F3

%%% Generate and plot time series
subplot(2,4,2) 
hold on
beta = 0.1;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); set(gca,'XTick',[0 2 4 6])
ylim([10 24]); set(gca,'YTick',[10 17 24])
xlabel({' '}) 
ylabel({' '}) 
title('b: F_3') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel C: F4

%%% Generate and plot time series
subplot(2,4,3) 
hold on
beta = 0.1;
eta = 0.5;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + eta*E - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); set(gca,'XTick',[0 2 4 6])
ylim([10 24]); set(gca,'YTick',[10 17 24])
xlabel({' '}) 
ylabel({' '}) 
title('c: F_4') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel D: F5

%%% Generate and plot time series
subplot(2,4,5) 
hold on
beta = 2;
eta = 0.5;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C + eta*E - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); set(gca,'XTick',[0 2 4 6])
ylim([0 32]); set(gca,'YTick',[0 16 32])
xlabel({'Time'}) 
ylabel({'Ecosystem function'}) 
title('d: F_5') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel E: F6

%%% Generate and plot time series
subplot(2,4,6) 
hold on
beta = 0.1;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C./(1 + gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); set(gca,'XTick',[0 2 4 6])
ylim([19 21]); set(gca,'YTick',[19 20 21])
xlabel({'Time'}) 
ylabel({' '}) 
title('e: F_6') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel F: F6

%%% Generate and plot time series
subplot(2,4,7) 
hold on
beta = 0.1;
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = (alpha + beta*C).*exp(-gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 9]); set(gca,'XTick',[0 3 6 9])
ylim([0 24]); set(gca,'YTick',[0 12 24])
xlabel({'Time'}) 
ylabel({' '}) 
title('f: F_7') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
