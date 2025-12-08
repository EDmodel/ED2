function plot_soilcarbon(dns,scp_t,scp_c,titlestr,soilcarb_pref,visible)

global fasz;
global scpname;
global scpcolor;
global nscp;

fasz_l = fasz+1;
fasz_c = fasz-0;

figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches');
set(gcf,'Position',[0.25 0.25 5 6]);

% Determine which pfts are part of the game


% Find total Soil Carbon
tot_scp_t = sum(scp_t,2);
tot_scp_c = sum(scp_c,2);


xticks = linspace(min(dns),max(dns),5);
xvec = datevec(xticks);
xticklabs = num2str(xvec(:,1));
xticklabs(end,:) = '    ';

ymax = 1.04 .* max([tot_scp_t;tot_scp_c]);

ax1 = axes('Position',[0.15 0.30 0.80 0.65],'FontSize',fasz_l);
hold on;
p_c=plot(dns,[scp_c(:,:) tot_scp_c],'LineStyle','--');
p_t=plot(dns,[scp_t(:,:) tot_scp_t],'LineStyle','-' );
for ip=1:(nscp+1)
   if (ip == (nscp+1))
      set(p_c(ip),'Color','k','LineWidth',1.75);
      set(p_t(ip),'Color','k','LineWidth',1.75);
   else
      set(p_c(ip),'Color',scpcolor(ip,:),'LineWidth',1.75);
      set(p_t(ip),'Color',scpcolor(ip,:),'LineWidth',1.75);
   end
end
hold off;
datetick;
grid on;
box on;
ylabel('Soil Carbon [kgC/m^2]','FontSize',fasz_l)
scpleg=[scpname 'Total'];
lhan=legend(ax1,scpleg,'Location','SouthOutside');
set(lhan,'Position',[0.35 0.05 0.35 0.15],'FontSize',fasz_c)
xlim([min(dns) max(dns)]);
ylim([0 ymax]);
set(gca,'XTick',xticks,'XTickLabel',xticklabs,'FontSize',fasz_l);
title('Soil Carbon (dashed is MainLine)','FontSize',fasz_l);


oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches','PaperPosition',newpos)
print('-depsc',sprintf('%s.eps',soilcarb_pref),'-r300');
print('-dpng' ,sprintf('%s.png',soilcarb_pref),'-r300');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)

%print('-dpng','-r100',pftsucc_img);
