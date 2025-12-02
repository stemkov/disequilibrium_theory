% AUTHOR:   Michael Cortez, Florida State University
% DATE:     December 2025
% PURPOSE:  Plot time series for formulas F3 & F4, assuming linearly increasing climate

% INPUTS:   Parameters for ecosystem function formula (alpha, beta, gamma, kappa)
%           Delta0 = Initial disequilibirum when climate change begins
%           C0 = initial climate
%           r = rate of change of climate
%           lambda = response rate of community

clear all

figure('Position',[10 10 1200 600]);
textsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel A

%%% Define parameters
tvals = 0:0.01:20;
alpha = 20;
beta = 1;
kappa = 0.5; 
gamma =4;
C0 = 2; 
r = 1;
lambda = 2;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 1 r*lambda 4 10]
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,4,1) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 10 20])
ylim([-20 30])
xlabel('time') 
ylabel({'Linear dependence'; 'form'; 'Ecosystem function (F_3)'}) 
title('a') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('\Delta_0>0, Q<0','\Delta_0>0, Q<0','\Delta_0>0, Q>0',...
    '\Delta_0=0, Q>0','\Delta_0<0, Q>0','\Delta_0<0, Q>0',...
    'position',[0.35 0.1 1 1],'FontSize', 13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel D

%%% Define parameters
% Uses same parameters as Panel A

%%% Generate and plot time series
subplot(2,4,5) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha +  beta*(1-exp(-C*kappa)) - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 10 20])
ylim([-20 30])
xlabel('time') 
ylabel({'Nonlinear dependence'; 'form'; 'Ecosystem function (F_3)'}) 
title('d') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('\Delta_0>0, Q<0','\Delta_0>0, Q<0','\Delta_0>0, Q>0',...
    '\Delta_0=0, Q>0','\Delta_0<0, Q>0','\Delta_0<0, Q>0',...
    'position',[0.35 0.1 1 1],'FontSize', 13)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel B

%%% Define parameters
tvals = 0:0.01:10;
alpha = 20;
eta = 4;
gamma =0.1;
C0 = 2; 
r = 1;
lambda = 2;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 1 r*lambda 4 10]
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,4,2) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + eta*E - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([-20 50])
xlabel('time') 
ylabel('Ecosystem function (F_4)') 
title('b') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel C

%%% Define parameters 
tvals = 0:0.01:10;
alpha = 20;
eta = 1;
gamma =4;
C0 = 2; 
r = 1;
lambda = 2;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 1 r*lambda 4 10]
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,4,3) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + eta*E - gamma.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([-20 50])
xlabel('time') 
ylabel('Ecosystem function (F_4)') 
title('c') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
