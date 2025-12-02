% AUTHOR:   Michael Cortez, Florida State University
% DATE:     December 2025
% PURPOSE:  Plot time series for formulas F6 & F7, assuming linearly increasing climate

% INPUTS:   Parameters for ecosystem function formula (alpha, beta, gamma, kappa)
%           Delta0 = Initial disequilibirum when climate change begins
%           C0 = initial climate
%           r = rate of change of climate
%           lambda = response rate of community

clear all

figure('Position',[10 10 900 600]);
textsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel A

%%% Define parameters
tvals = 0:0.01:10;
alpha = 1;
beta = 1;
gamma =1;
kappa = .5; 
C0 = 4; 
r = 1;
lambda = 1;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -2 0.7 r*lambda 2 10]
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,3,1) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C./(1 + gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([0 10])
xlabel('time') 
ylabel({'Linear dependence'; 'form'; 'Ecosystem function (F_6)'}) 
title('a') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('\Delta_0>0, Q<0','\Delta_0>0, Q<0','\Delta_0>0, Q>0',...
    '\Delta_0=0, Q>0','\Delta_0<0, Q>0','\Delta_0<0, Q>0',...
    'position',[0.35 0.1 1 1],'FontSize', 13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel C

%%% Define parameters
% Uses same parameters as Panel A

%%% Generate and plot time series
subplot(2,3,4) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*(1-exp(-C*kappa))./(1 + gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([0 3])
xlabel('time') 
ylabel({'Nonlinear dependence'; 'form'; 'Ecosystem function (F_6)'}) 
title('c') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel B

%%% Define parameters 
tvals = 0:0.01:10;
alpha = 1;
beta = 1;
gamma =1;
C0 = 4; 
r = 1;
lambda = 1;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -2 0.7 r*lambda 2 10]
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,3,2) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = (alpha + beta*C).*exp(-gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([0 10])
xlabel('time') 
ylabel('Ecosystem function (F_7)') 
title('b') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('\Delta_0>0, Q<0','\Delta_0>0, Q<0','\Delta_0>0, Q>0',...
    '\Delta_0=0, Q>0','\Delta_0<0, Q>0','\Delta_0<0, Q>0',...
    'position',[0.35 0.1 1 1],'FontSize', 13)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel D

%%% Define parameters
% Uses same parameters as Panel B

%%% Generate and plot time series
subplot(2,3,5) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = (alpha + beta*(1-exp(-C*kappa))).*exp(-gamma.*(E-C).^2);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 5 10])
ylim([0 3])
xlabel('time') 
ylabel('Ecosystem function (F_7)') 
title('d') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
