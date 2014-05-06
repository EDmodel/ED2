function plot_fluxes(flux_tp,flux_cp,flux_td,flux_cd,flux_tm,flux_cm,...
    varnames,varunits,phrs,titlestr,fluxes_img,visible)

global fasz;

ntmp = 5;

fasz_l=fasz+1;

figure('visible',visible);
set(gcf,'PaperPositionMode','manual',...
    'Units','inches','Position',[0.25 0.25 10.5 6]);


dx = 0.21;bx = 0.07;mx = 0.05;
dy = 0.24;by = 0.08;my = 0.05;
dx = ((0.97-bx)./ntmp)-mx;
       

for ip=1:ntmp
    
    pos=ntmp-ip;
  
    ax1 = axes;
    set(ax1,'Position',[bx+pos*dx+pos*mx by dx dy],'FontSize',fasz_l);

    minxy = min([min(flux_tp(:,ip)) min(flux_cp(:,ip))]);
    maxxy = max([max(flux_tp(:,ip)) max(flux_cp(:,ip))]);
    plot(sort(flux_tp(:,ip)),sort(flux_cp(:,ip)),'o');
    xlim([minxy maxxy]);
    ylim([minxy maxxy]);
    grid on;
    box on;
    
    xlabel('');
    if(ip==ntmp);
        ylabel(sprintf('Control/Test\nQuantiles'));
    else
        ylabel('');
    end
    
    ax2 = axes;
    set(ax2,'Position',[bx+pos*dx+pos*mx by+my+dy dx dy],'FontSize',fasz_l);
    hold on;
    plot(phrs,flux_td(:,ip),'o','Color','b');
    plot(phrs,flux_cd(:,ip),'+','Color','r');
    
    miny = min([min(flux_td(:,ip)) min(flux_cd(:,ip))]);
    maxy = max([max(flux_td(:,ip)) max(flux_cd(:,ip))]);
    
    miny=miny-0.03*abs(maxy);
    maxy=maxy+0.03*abs(maxy);
    
    ylim([miny maxy]);
    xlim([1 24]);
    set(gca,'XTick',[6 12 18]);
    hold off;
    grid on;
    box on;
    %if(ip==1);title(titlestr,'FontSize',fasz);end;
    if(ip==ntmp);ylabel('Hourly Means','FontSize',fasz_l);end;

    ax3 = axes;
    set(ax3,'Position',[bx+pos*dx+pos*mx by+2*my+2*dy dx dy],...
        'FontSize',fasz_l);
    hold on;
    plot(1:12,flux_tm(:,ip),'o','Color','b');
    plot(1:12,flux_cm(:,ip),'+','Color','r');
    hold off;
    xlim([1 12]);
    set(gca,'XTick',[3 6 9],'XTicklabel',{'Mar','Jun','Sep'},'FontSize',fasz_l);
    grid on;
    box on;
    if(ip==ntmp);ylabel('Monthly Means','FontSize',fasz_l);end;
    title(sprintf('%s\n%s',varnames{ip},varunits{ip}),'FontSize',fasz_l);
end

ax4=axes;
set(ax4,'Position',[0.0 0.02 1.0 0.05]);
axis off;
text(0.5,0.05,'Test: o       Main:+','HorizontalAlignment','center','FontSize',fasz_l);


oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc',fluxes_img,'-r200');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)


%print('-dpng','-r300',fluxes_img);
%print('-djpeg','-r300','flux.jpg');

