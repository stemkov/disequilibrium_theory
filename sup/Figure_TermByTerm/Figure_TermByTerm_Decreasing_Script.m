% AUTHOR:   Michael Cortez, Florida State University
% DATE:     December 2025
% PURPOSE:  Plot time series for each term, assuming linearly decreasing climate

% INPUTS:   Parameters for ecosystem function formula (alpha, beta, gamma, kappa)
%           Delta0 = Initial disequilibirum when climate change begins
%           C0 = initial climate
%           r = rate of change of climate
%           lambda = response rate of community


clear all

figure('Position',[10 10 1200 900]);
textsize = 14;

%% Common Parameters - Rows 1 and 2
tvals = 0:0.01:10;
r = 1;
C0 = 5;
C = C0 - r*tvals;
Colors = {'b','c', 'k', 'm','r'};
lambda = 1;
Delta0 = [3, 1, 0, -0.5, -2];
F = zeros(length(Delta0),length(tvals));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel A: E
subplot(3,5,1) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = E;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([0 8]); set(gca,'YTick',[0 2 4 6 8])
xlabel('') 
ylabel('E') 
title('a: E') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')
legend('$\Delta_0>0$, $\bar{Q}<0$','$\Delta_0>0$, $\bar{Q}<0$',...
    '$\Delta_0$=0, $\bar{Q}<0$','$\Delta_0<0$, $\bar{Q}<0$','$\Delta_0<0$, $\bar{Q}>0$',...
    'position',[0.25 0.305 1 1],'FontSize', 16, 'Interpreter', 'Latex',...
    'FontWeight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel B: -|E-C|
subplot(3,5,2) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -abs(E-C);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-3 0]); set(gca,'YTick',[-3 -2 -1 0])
ylabel('-|\Delta|') 
xlabel('') 
title('b: -|\Delta|') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel C: -(E-C)^2
subplot(3,5,3) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-4 0]); set(gca,'YTick',[-4 -3 -2 -1 0])
ylabel('-\Delta^2') 
xlabel('') 
title('c: -\Delta^2') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel D: EC
subplot(3,5,6) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = C.*E;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 10]); 
set(gca,'XTick',[0 2 4 6 8 10])
ylim([0 40]); set(gca,'YTick',[0 10 20 30 40])
ylabel('CE') 
xlabel('time') 
title('d: EC') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel E: -C|E-C|
subplot(3,5,7) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -C.*abs(E-C);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-8 0]); set(gca,'YTick',[-10 -8 -6 -4 -2 0])
ylabel('-C|\Delta|') 
xlabel('time') 
title('e: -C|\Delta|') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel F: C(E-C)^2
subplot(3,5,8) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -C.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-8 0]); set(gca,'YTick',[-10 -8 -6 -4 -2 0])
ylabel('-C\Delta^2') 
xlabel('time') 
title('f: -C\Delta^2') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel G: C|E-C|
subplot(3,5,9) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -E.*abs(E-C);
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-8 0]); set(gca,'YTick',[-10 -8 -6 -4 -2 0])
ylabel('-E|\Delta|') 
xlabel('time') 
title('g: -E|\Delta|') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Panel H: C(E-C)^2
subplot(3,5,10) 
hold on
for ii = 1:length(Delta0)
    %%% Compute time series for community climate equilibrium   
    E = C + r*lambda + (-r*lambda - Delta0(ii))*exp(-tvals/lambda);
    %%% Compute time series for ecosystem function
    F(ii,:) = -E.*(E-C).^2;
    plot(tvals,F(ii,:),'-','LineWidth',1.5,'color',Colors{ii})
end
xlim([0 6]); 
set(gca,'XTick',[0 2 4 6])
ylim([-8 0]); set(gca,'YTick',[-10 -8 -6 -4 -2 0])
ylabel('-E\Delta^2') 
xlabel('time') 
title('h: -E\Delta^2') 
set(gca,'Box','on', 'Layer','top','LineWidth', 1)
set(gca,'ycolor','k','FontSize',textsize,'FontWeight','Bold')


