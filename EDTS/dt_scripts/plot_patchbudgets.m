function plot_patchbudgets(dnv,cbuds_t,ebuds_t,wbuds_t,...
                                cstor_t,estor_t,wstor_t,...
                                cbuds_c,ebuds_c,wbuds_c,...
                                cstor_c,estor_c,wstor_c,...
                                cbudg_outfile, ...
                                ebudg_outfile, ...
                                wbudg_outfile, ...
                                visible);
global fasz;



cbud_cols=[203     0   61;     ... %Residual
           127   127  127;     ... %Delta S
           151   228  125;     ... %Density
            16   124  146;     ... %Veg Dyn
           248   120   86;     ... %Net Seed
            27   162  247]./255.;  %Eddy flux

ebud_cols=[203     0   61;    ... % Residual
           127   127  127;    ... % Delta S
            43   210  219;    ... % Rainfall
           248   120   86;    ... % Net Rad
           151   228  125;    ... % Density
           204   202   61;    ... % Pressure
            16   124  146;    ... % Veg Dyn
            27   162  247;    ... % Eddy flux
           129    31  158]./255.; % Runoff


wbud_cols=[203     0   61;    ... % Residual
           127   127  127;    ... % Delta S
            43   210  219;    ... % Rainfall
           151   228  125;    ... % Density
            27   162  247;    ... % Eddy flux
            16   124  146;    ... % Veg Dyn
           129    31  158]./255.; % Runoff

cbud_names = {'Residual','\Delta S','Density','Veg Dyn','Net Seed','Eddy Flux'};
ebud_names = {'Residual','\Delta S','Rainfall','R_{net}','Density','Pressure','Veg Dyn','Eddy Flux','Runoff'};
wbud_names = {'Residual','\Delta S','Rainfall','Density','Veg Dyn','Eddy Flux','Runoff'};

error_cols=[230    92   23; ...    % Test
             41   150  204]./255.; % Test

error_names={'Test','Main'};

fasz_l = fasz+1;
fasz_c = fasz-2;

nebud = size(ebud_cols,1);
nwbud = size(wbud_cols,1);
ncbud = size(cbud_cols,1);


% Energy balance plots for patch 1

by = 0.08; my = 0.03;
bx = 0.10; mx = 0.04;
dy = 0.26;
dx = 0.8;

[~,~,npatch] = size(ebuds_t);

npatch = min(size(nebud,1)-1,npatch);
dtfac  = (dnv(2)-dnv(1))*86400.0; %Time integrator (seconds)

% =========================================================================
% Energy Budget
% =========================================================================
figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches','PaperSize',[11 8.5]);
set(gcf,'Position',[0.25 0.25 10 7.5]);

ax1 = axes('Position',[bx by+2*(dy+my) dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,1e-6*dtfac*cumsum(ebuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,1e-6*dtfac*cumsum(ebuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:nebud
    set(p_c(ip),'Color',ebud_cols(ip,:));
    set(p_t(ip),'Color',ebud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Energy Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[GJ/m2]'),'FontSize',fasz_l);
legend(ebud_names,'Location','NorthWest','FontSize',fasz_c);
set(ax1,'XtickLabel',{});

ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,ebuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,ebuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:nebud
    set(p_c(ip),'Color',ebud_cols(ip,:));
    set(p_t(ip),'Color',ebud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[W/m2]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});

ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p_t=plot(dnv,1e-6*dtfac*cumsum(ebuds_t(:,1,1),1),...
    'LineStyle','-','LineWidth',1.5);
p_c=plot(dnv,1e-6*dtfac*cumsum(ebuds_c(:,1,1),1),...
    'LineStyle','--','LineWidth',1.5);
set(p_t(1),'Color',error_cols(1,:));
set(p_c(1),'Color',error_cols(2,:));
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[GJ/m2]'),'FontSize',fasz_l);
legend(error_names,'Location','NorthWest','FontSize',fasz_c);

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc', ebudg_outfile, '-r200');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)

% =========================================================================
% Water Budget
% =========================================================================

idaysec = 1/86400;

figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches','PaperSize',[11 8.5]);
set(gcf,'Position',[0.25 0.25 10 7.5]);

ax1 = axes('Position',[bx by+2*(dy+my) dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,idaysec*dtfac*cumsum(wbuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,idaysec*dtfac*cumsum(wbuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:nwbud
    set(p_c(ip),'Color',wbud_cols(ip,:));
    set(p_t(ip),'Color',wbud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Water Mass Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[kg/m2]'),'FontSize',fasz_l);
legend(wbud_names,'Location','NorthWest','FontSize',fasz_c);
set(ax1,'XtickLabel',{});

ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,3600*idaysec*wbuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,3600*idaysec*wbuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:nwbud
    set(p_c(ip),'Color',wbud_cols(ip,:));
    set(p_t(ip),'Color',wbud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[kg/m2/hr]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});


ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p_t=plot(dnv,3600*idaysec*cumsum(wbuds_t(:,1,1),1),...
    'LineStyle','-','LineWidth',1.5);
p_c=plot(dnv,3600*idaysec*cumsum(wbuds_c(:,1,1),1),...
    'LineStyle','--','LineWidth',1.5);
set(p_t(1),'Color',error_cols(1,:));
set(p_c(1),'Color',error_cols(2,:));
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[kg/m2]'),'FontSize',fasz_l);
legend(error_names,'Location','NorthWest','FontSize',fasz_c);

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc', wbudg_outfile, '-r200');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)

% =========================================================================
% Carbon Budget
% =========================================================================
idaysec = 1/86400;

figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches','PaperSize',[11 8.5]);
set(gcf,'Position',[0.25 0.25 10 7.5]);

ax1 = axes('Position',[bx by+2*(dy+my) dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,dtfac*cumsum(cbuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,dtfac*cumsum(cbuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:ncbud
    set(p_c(ip),'Color',cbud_cols(ip,:));
    set(p_t(ip),'Color',cbud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Carbon Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[kgC/m2]'),'FontSize',fasz_l);
legend(cbud_names,'Location','NorthWest','FontSize',fasz_c);
set(ax1,'XtickLabel',{});


ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p_c=plot(dnv,cbuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
p_t=plot(dnv,cbuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
for ip=1:ncbud
    set(p_c(ip),'Color',1000.*idaysec*cbud_cols(ip,:));
    set(p_t(ip),'Color',1000.*idaysec*cbud_cols(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[gC/m2/dy]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});


ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p_t=plot(dnv,dtfac*cumsum(cbuds_t(:,1,1),1),...
    'LineStyle','-','LineWidth',1.5);
p_c=plot(dnv,dtfac*cumsum(cbuds_c(:,1,1),1),...
    'LineStyle','--','LineWidth',1.5);
set(p_t(1),'Color',error_cols(1,:));
set(p_c(1),'Color',error_cols(2,:));
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[kgC/m2]'),'FontSize',fasz_l);
legend(error_names,'Location','NorthWest','FontSize',fasz_c);



oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc', cbudg_outfile, '-r200');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)
