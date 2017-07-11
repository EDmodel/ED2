function plot_patchbudgets(dnv,cbuds_t,ebuds_t,wbuds_t,...
                                cstor_t,estor_t,wstor_t,...
                                cbuds_c,ebuds_c,wbuds_c,...
                                cstor_c,estor_c,wstor_c,...
                                cbudg_outfile, ...
                                ebudg_outfile, ...
                                wbudg_outfile, ...
                                visible);
global fasz;
global pftcolor;


fasz_l = fasz+1;


% Energy balance plots for patch 1

by = 0.08; my = 0.03;
bx = 0.10; mx = 0.04;
dy = 0.26;
dx = 0.8;

ebud_names = {'Residual','\Delta S','Precip','R_{net}','Density','Pressure','Atmosphere','Drainage','Runoff'};
cbud_names = {'Residual','\Delta S','Nep','Density','Atsmosphere'};
wbud_names = {'Residual','\Delta S','Precip','Density','Atmosphere','Drainage','Runoff'};

[~,~,npatch] = size(ebuds_t);

npatch = min(size(pftcolor,1)-1,npatch);

for ip=1:npatch
    patch_names{ip} = sprintf('p%i',ip);
end

dtfac = (dnv(2)-dnv(1))*86400.0; %Time integrator (seconds)

% =========================================================================
% Energy Budget
% =========================================================================
figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches','PaperSize',[11 8.5]);
set(gcf,'Position',[0.25 0.25 10 7.5]);

ax1 = axes('Position',[bx by+2*(dy+my) dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,1e-6*dtfac*cumsum(ebuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,1e-6*dtfac*cumsum(ebuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:9
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Energy Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[GJ/m2]'),'FontSize',fasz_l);
legend(ebud_names,'Location','NorthWest','FontSize',fasz_l);
set(ax1,'XtickLabel',{});

ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,ebuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,ebuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:9
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[W/m2]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});

ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,ebuds_t(:,1,:),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,ebuds_c(:,1,:),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:npatch
    set(p1(ip),'Color',pftcolor(ip+1,:));
    set(p2(ip),'Color',pftcolor(ip+1,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[W/m2]'),'FontSize',fasz_l);
legend(patch_names,'Location','NorthWest','FontSize',fasz_l);

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
p1=plot(dnv,idaysec*dtfac*cumsum(wbuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,idaysec*dtfac*cumsum(wbuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:7
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Water Mass Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[kg/m2]'),'FontSize',fasz_l);
legend(wbud_names,'Location','NorthWest','FontSize',fasz_l);
set(ax1,'XtickLabel',{});

ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,3600*idaysec*wbuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,3600*idaysec*wbuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:7
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[kg/m2/hr]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});


ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,3600*idaysec*wbuds_t(:,1,:),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,3600*idaysec*wbuds_c(:,1,:),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:npatch
    set(p1(ip),'Color',pftcolor(ip+1,:));
    set(p2(ip),'Color',pftcolor(ip+1,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[kg/m2/hr]'),'FontSize',fasz_l);
legend(patch_names,'Location','NorthWest','FontSize',fasz_l);

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

figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches','PaperSize',[11 8.5]);
set(gcf,'Position',[0.25 0.25 10 7.5]);

ax1 = axes('Position',[bx by+2*(dy+my) dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,dtfac*cumsum(cbuds_t(:,:,1),1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,dtfac*cumsum(cbuds_c(:,:,1),1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:5
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
title(sprintf('Carbon Budget Terms\n(dashed lines for mainline)'),...
    'Fontsize',fasz_l);
ylabel(sprintf('Integrated\n[umol/m2]'),'FontSize',fasz_l);
legend(cbud_names,'Location','NorthWest','FontSize',fasz_l);
set(ax1,'XtickLabel',{});


ax2 = axes('Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,cbuds_t(:,:,1),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,cbuds_c(:,:,1),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:5
    set(p1(ip),'Color',pftcolor(ip,:));
    set(p2(ip),'Color',pftcolor(ip,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Instananeous\n[umol/m2/s]'),'FontSize',fasz_l);
set(ax2,'XtickLabel',{});


ax3 = axes('Position',[bx by dx dy],'FontSize',fasz_l);
hold on;
p1=plot(dnv,cbuds_t(:,1,:),...
    'LineStyle','-','LineWidth',1.5);
p2=plot(dnv,cbuds_c(:,1,:),...
    'LineStyle','--','LineWidth',1.5);
for ip=1:npatch
    set(p1(ip),'Color',pftcolor(ip+1,:));
    set(p2(ip),'Color',pftcolor(ip+1,:));
end
hold off;
datetick;
grid on;
box on;
ylabel(sprintf('Error \n[umol/m2/s]'),'FontSize',fasz_l);
legend(patch_names,'Location','NorthWest','FontSize',fasz_l);



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
