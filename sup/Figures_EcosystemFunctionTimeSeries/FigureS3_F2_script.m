% AUTHOR:   Michael Cortez, Florida State University
% DATE:     December 2025
% PURPOSE:  Plot time series for formula F2, assuming linearly increasing climate

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
tvals = 0:0.01:8;
alpha = 10;
beta = 1;
gamma =1;
kappa = 0.5;
C0 = 1; 
r = 1;
lambda = 0.8;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 0 r*lambda+0.2 2 10]+0.2
F = zeros(length(Delta0),length(tvals));

%%% Generate and plot time series
subplot(2,4,1) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([-10 15])
xlabel('time') 
ylabel({'Linear dependence'; 'form'; 'Ecosystem function (F_2)'}) 
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
    F(ii,:) = alpha + beta*(1-exp(-C*kappa)) - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([-10 15])
xlabel('time') 
ylabel({'Nonlinear dependence'; 'form'; 'Ecosystem function (F_2)'}) 
title('d') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel B

%%% Define parameters
tvals = 0:0.01:8;
alpha = 12;
beta = 1;
gamma =1;
C0 = 0; 
r = 1;
lambda = 0.8;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 0 r*lambda+0.2 2 10]+0.2
F = zeros(length(Delta0),length(tvals));

subplot(2,4,2) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([0 15])
xlabel('time') 
ylabel('') 
title('b') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel E

%%% Define parameters
% Uses same parameters as Panel B

%%% Generate and plot time series
subplot(2,4,6)
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*(1-exp(-C*kappa)) - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([0 15])
xlabel('time') 
ylabel('') 
title('e') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel C

%%% Define parameters
tvals = 0:0.01:8;
alpha = 20;
beta = 1;
gamma =1.5;
C0 = 0; 
r = 1;
lambda = 1;
C = C0 + r*tvals;
Colors = {'b',"#80B3FF",'c', 'k', 'm','r'};
Delta0 = r*lambda - [-5 -1 0 r*lambda 2 10]
F = zeros(length(Delta0),length(tvals));


subplot(2,4,3) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*C - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([0 30])
xlabel('time') 
ylabel('') 
title('c') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel F

%%% Define parameters
% Uses same parameters as Panel D

%%% Generate and plot time series
subplot(2,4,7)
hold on

for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium    
    E = C - r*lambda + (r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = alpha + beta*(1-exp(-C*kappa)) - gamma.*C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 max(tvals)]); set(gca,'XTick',[0 4 8])
ylim([0 30])
xlabel('time') 
ylabel('') 
title('f') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
