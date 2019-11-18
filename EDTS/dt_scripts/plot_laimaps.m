function plot_laimaps(usepft,lai_gt,lai_gc,lon_gc,...
    lat_gc,npoly,laimap_pref,visible,grid_name)

global fasz;
global pftname;
global actcmap;
global diffcmap;
load rywcbmap.mat;
load wygmap.mat;

[lonpcrns,latpcrns] = approx_patch_4corners(double(lat_gc),double(lon_gc));

% Match the validation grid to the model grid

nupft = sum(usepft);

fid=fopen('americadosul_unf');
C=textscan(fid,'%n%n');
fclose(fid);

ids=find(C{1}==999);
nbnds = length(ids);

prev=1;
for i=1:nbnds
next=ids(i)-1;
geodata(i) = struct('lon',C{1}(prev:next),'lat',C{2}(prev:next)); 
prev=ids(i)+1;
end


figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches');
set(gcf,'Position',[0.25 0.25 min((nupft+1)*3,10.0) 6]);

by = 0.15; my = 0.1;
bx = 0.08; mx = 0.01;
dy = 0.3;
dx = ((0.97-bx)./(nupft+1))-mx;

minlon = min(lon_gc);
maxlon = max(lon_gc);
minlat = min(lat_gc);
maxlat = max(lat_gc);
dlon = maxlon-minlon;
dlat = maxlat-minlat;
minlon = minlon - 0.15*dlon;
maxlon = maxlon + 0.15*dlon;
minlat = minlat - 0.15*dlat;
maxlat = maxlat + 0.15*dlat;

% Title

ax = axes;
set(ax,'Position',[0 0.95 1 0.05]);
axis off;
text(0.5,0.45,sprintf('LAI - %s',grid_name),...
    'FontSize',fasz,'HorizontalAlignment','center');


% LEFT - Total LAI

patch_lai_dgt = zeros(4,npoly);
patch_lai_gc = zeros(4,npoly);
for ipy=1:npoly
    tot_lai_gc = sum(lai_gc(ipy,:));
    tot_lai_gt = sum(lai_gt(ipy,:));
    if (tot_lai_gc == 0. && tot_lai_gt == 0.)
       patch_lai_dgt(:,ipy) = 0.0;
    else
       patch_lai_dgt(:,ipy) = 200.*(tot_lai_gt-tot_lai_gc)...
                           ./ (tot_lai_gc + tot_lai_gt);
    end
    patch_lai_gc(:,ipy) = sum(lai_gc(ipy,:));
end
minc = 0.0;
maxc = max(max(patch_lai_gc));

ax1 = axes;
set(ax1,'Position',[bx by+dy+my dx dy],'FontSize',fasz);
hold on;
patch(lonpcrns,latpcrns,patch_lai_gc);
colormap(actcmap);
grid on; box on;
caxis([minc maxc]);
shading flat;
ylabel('Main','Fontsize',12);
set(gca,'XtickLabel',{});
cobar=colorbar('South','Position',[bx+0.05*dx by+dy+0.02 0.9*dx 0.03],'FontSize',fasz);
set(cobar,'XTick',[minc,maxc]);
for b=1:nbnds
plot(geodata(b).lon,geodata(b).lat,'Color',[0.9 0.5 0.5],'LineWidth',1.0);
end
hold off;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
title('Total LAI [m2/m2]','FontSize',fasz);
freezeColors;
cbfreeze;


maxdc = max([1,max(max(abs(patch_lai_dgt)))]);
mindc = -maxdc;


ax2 = axes;
set(ax2,'Position',[bx by dx dy],'FontSize',fasz);
hold on;
patch(lonpcrns,latpcrns,patch_lai_dgt);
colormap(diffcmap);
grid on; box on;
caxis([mindc maxdc]);
shading flat;
colorbar('South','Position',[bx+0.05*dx 0.05 0.9*dx 0.03],'FontSize',fasz)
for b=1:nbnds
plot(geodata(b).lon,geodata(b).lat,'Color',[0.9 0.5 0.5],'LineWidth',1.0);
end
hold off;
xlim([minlon maxlon])
ylim([minlat maxlat]);
ylabel('200(T-M)/(T+M)','FontSize',fasz);
freezeColors;
cbfreeze;

ipfts=find(usepft>0);

% Partitions


for ip=1:numel(ipfts)

ipft=ipfts(ip);    

patch_lai_dgt = zeros(4,npoly);
patch_lai_gc = zeros(4,npoly);

for ipy=1:npoly
    if (lai_gc(ipy,ipft) == 0. && lai_gt(ipy,ipft) == 0.)
       patch_lai_dgt(:,ipy) = 0.;
    else
       patch_lai_dgt(:,ipy) = 200. * (lai_gt(ipy,ipft)-lai_gc(ipy,ipft))...
                           ./ (lai_gt(ipy,ipft)+lai_gc(ipy,ipft));
    end
%    patch_lai_dgt(:,ipy) = lai_gt(ipy,ipft)-lai_gc(ipy,ipft);
    patch_lai_gc(:,ipy) = lai_gc(ipy,ipft);
end
minc = 0.0;
maxc = max(max(patch_lai_gc));


ax = axes; %#ok<LAXES>
set(ax,'Position',[bx+ip*(dx+mx) by+my+dy dx dy],'FontSize',fasz);
hold on;
patch(lonpcrns,latpcrns,patch_lai_gc);
colormap(actcmap);
grid on; box on;
caxis([minc maxc]);
shading flat;
title(sprintf('%s',pftname{ipft}),'Fontsize',fasz);
colorbar('South','Position',[bx+ip*(dx+mx)+0.05*dx by+dy+0.02 0.9*dx 0.03],'FontSize',fasz)
set(gca,'XtickLabel',{});
set(gca,'YtickLabel',{});
freezeColors;
cbfreeze;

for b=1:nbnds
plot(geodata(b).lon,geodata(b).lat,'Color',[0.9 0.5 0.5],'LineWidth',1.0);
end
hold off;
xlim([minlon maxlon]);
ylim([minlat maxlat]);

maxdc = max([1,max(max(abs(patch_lai_dgt)))]);
mindc = -maxdc;

%maxdc = max([0.001*maxc  ,max(max(abs(patch_lai_dgt)))]);
%mindc = -maxdc;

ax = axes; %#ok<LAXES>
set(ax,'Position',[bx+ip*(dx+mx) by dx dy],'FontSize',fasz);
hold on;
patch(lonpcrns,latpcrns,patch_lai_dgt);
colormap(diffcmap);
grid on; box on;
caxis([mindc maxdc]);
shading flat;
colorbar('South','Position',[bx+ip*(dx+mx)+0.05*dx 0.05 0.9*dx 0.03],'FontSize',fasz)
set(gca,'YtickLabel',{});
for b=1:nbnds
plot(geodata(b).lon,geodata(b).lat,'Color',[0.9 0.5 0.5],'LineWidth',1.0);
end
hold off;
xlim([minlon maxlon]);
ylim([minlat maxlat]);
freezeColors;
cbfreeze;

end


oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-depsc', sprintf('%s.eps',laimap_pref), '-r300');
print('-dpng', sprintf('%s.png',laimap_pref), '-r300');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)



