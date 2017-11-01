function plot_succession(dns,agb_t,agb_c,lai_t,lai_c,...
    titlestr,pftsucc_pref,visible)

global fasz;
global pftshort;
global pftcolor;
global npft;

fasz_l = fasz+1;

figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches');
set(gcf,'Position',[0.25 0.25 9 6.5]);

bx = 0.1; by=0.1;
mx=0.03;  my=0.05;
dx = 0.3; dy = 0.35;

% Determine which pfts are part of the game

clear upft;
ip=0;
for ipft=1:npft
    if(~isempty(find(agb_t(:,ipft)>0)) || ~isempty(find(agb_c(:,ipft)))) %#ok<*EFIND>
        ip=ip+1;
        upft(ip) = ipft; %#ok<AGROW>
    end
end

if(numel(upft)<1) %Force at least some plotting, C4
    upft(1)=1;
end


% Find total AGB and LAI
tot_agb_t = sum(agb_t,2);
tot_agb_c = sum(agb_c,2);
tot_lai_t = sum(lai_t,2);
tot_lai_c = sum(lai_c,2);

%-- AGB --

xticks = linspace(min(dns),max(dns),5);
xvec = datevec(xticks);
xticklabs = num2str(xvec(:,1));
xticklabs(end,:) = '    ';

ymax = 1.04 .* max([tot_agb_t;tot_agb_c]);

ax1 = axes;
set(ax1,'Position',[bx by+dy+my dx dy],'FontSize',fasz_l);
ph=plot(dns,[agb_t(:,upft) tot_agb_t]);
for ip=1:(numel(upft)+1)
  if (ip == (numel(upft)+1))
     set(ph(ip),'Color','k','LineWidth',1.75);
  else
     set(ph(ip),'Color',pftcolor(upft(ip),:),'LineWidth',1.25);
  end
end
datetick;
grid on;
box on;
ylabel('AGB [kgC/m^2]','FontSize',fasz_l)
set(gca,'XTick',xticks,'XTickLabel',{});
xlim([min(dns) max(dns)]);
ylim([0 ymax]);
title('Test','FontSize',fasz_l);

ax2 = axes;
set(ax2,'Position',[bx+dx+mx by+dy+my dx dy],'FontSize',fasz_l);
ph=plot(dns,[agb_c(:,upft) tot_agb_c]);
for ip=1:(numel(upft)+1)
  if (ip == (numel(upft)+1))
     set(ph(ip),'Color','k','LineWidth',1.75);
  else
     pftnow = upft(ip);
     set(ph(ip),'Color',pftcolor(upft(ip),:),'LineWidth',1.25);
  end
end
pftleg=[pftshort 'Total'];
lhan=legend(ax2,pftleg{[upft npft+1]},'Location','East');
set(lhan,'Position',[bx+2*dx+3*mx 0.35 0.1 0.3],'FontSize',fasz_l)
datetick;
grid on;
box on;
ylim([0 ymax]);
xlim([min(dns) max(dns)]);
set(gca,'XTick',xticks,'XTickLabel',{});
set(gca,'YTickLabel',{});
title('MainLine','FontSize',fasz_l);


ymax = 1.04 .* max([tot_lai_t;tot_lai_c]);

ax3 = axes;
set(ax3,'Position',[bx by dx dy],'FontSize',fasz_l);
ph=plot(dns,[lai_t(:,upft) tot_lai_t]);
for ip=1:(numel(upft)+1)
  if (ip == (numel(upft)+1))
     set(ph(ip),'Color','k','LineWidth',1.75);
  else
     pftnow = upft(ip);
     set(ph(ip),'Color',pftcolor(upft(ip),:),'LineWidth',1.25);
  end
end
datetick;
grid on;
box on;
ylabel('LAI [m^2/m^2]')
xlim([min(dns) max(dns)]);
set(gca,'XTick',xticks,'XTickLabel',xticklabs,'FontSize',fasz_l);
ylim([0 ymax]);

ax4 = axes;
set(ax4,'Position',[bx+dx+mx by dx dy],'FontSize',fasz_l);
ph=plot(dns,[lai_c(:,upft) tot_lai_c]);
for ip=1:(numel(upft)+1)
  if (ip == (numel(upft)+1))
     set(ph(ip),'Color','k','LineWidth',1.75);
  else
     pftnow = upft(ip);
     set(ph(ip),'Color',pftcolor(upft(ip),:),'LineWidth',1.25);
  end
end
datetick;
grid on;
box on;
xlim([min(dns) max(dns)]);
set(gca,'XTick',xticks,'XTickLabel',xticklabs,'FontSize',fasz_l);
set(gca,'YTickLabel',{});
ylim([0 ymax]);

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc',sprintf('%s.eps',pftsucc_pref),'-r300');
print('-dpng' ,sprintf('%s.png',pftsucc_pref),'-r300');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)

%print('-dpng','-r100',pftsucc_img);
