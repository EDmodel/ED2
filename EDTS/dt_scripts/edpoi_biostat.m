function edpoi_biostat(filename,ipy,poistr,figname,visible,type)

global fasz;
global pftname;
global pftshort;

% ========================================================
%
% Provide an ed state file. Define which polygon to use.
% Also, if this is part of a set, give it a name poistr.


ndbh = 11;
npft = 17;

% Get the vertical canopy layering
zmin = 0.40;
zmax = 1.0;
hmax = 36.0;
stretch=1.1;
[zzbot zztop] = canopy_layers(zmin,zmax,hmax,stretch);

ymax = hmax+2;
xmax = 0.5;
greeniness = [0.75  0.9   0.2;  %C4            (1)
	      0.2   0.8   0.2;
	      0.2   0.6   0.2;
	      0.2   0.4   0.2;
	      0.7   0.8   0.7;  %C3 Temperate  (5)
	      0.8   0.9   0  ;  %(6)
	      0.7   0.9   0  ;  %(7)
	      0.6   0.9   0  ;  %(8)
	      0.5   0.9   0  ;  %(9)
	      0.4   0.9   0  ;  %(10)
	      0.5   0.7   0  ;  %(11)
	      0.4   0.7   0  ;  %(12)
	      0.3   0.7   0  ;  %(13)
	      0.2   0.7   0  ;  %(14)
	      0.1   0.7   0  ;  %(15)
	      0.1   0.99  0.8;  %(16)
	      0.5   0.7   0.5;  %(17)
	     ];

% Read in raw data
% ========================================================
paco_id     = hdf5read(filename,'/PACO_ID');
sipa_id     = hdf5read(filename,'/SIPA_ID');
pysi_id     = hdf5read(filename,'/PYSI_ID');
paco_n      = hdf5read(filename,'/PACO_N');
sipa_n      = hdf5read(filename,'/SIPA_N');
pysi_n      = hdf5read(filename,'/PYSI_N');
npoly       = hdf5read(filename,'/NPOLYGONS_GLOBAL');
lat         = hdf5read(filename,'/LATITUDE');
lon         = hdf5read(filename,'/LONGITUDE');
if(type==1)
  agb_raw     = hdf5read(filename,'/AGB_PY');
  ba_raw      = hdf5read(filename,'/BASAL_AREA_PY');
 else
   ba_raw      = hdf5read(filename,'/BASAL_AREA_PY');
   agb_raw     = hdf5read(filename,'/AGB_PY');
end

area_si     = hdf5read(filename,'/AREA_SI');
area_pa     = hdf5read(filename,'/AREA');
age_pa      = hdf5read(filename,'/AGE');
slz         = hdf5read(filename,'/SLZ');
dist_pa     = hdf5read(filename,'/DIST_TYPE');

if(sum(paco_n)>0)
agb_co      = hdf5read(filename,'/AGB_CO');
height      = hdf5read(filename,'/HITE');
pft_co      = hdf5read(filename,'/PFT');
lai_co      = hdf5read(filename,'/LAI_CO');
nplant      = hdf5read(filename,'/NPLANT');
end



agb_pa = zeros(sum(sipa_n),1);

% Process data

lai_pft = zeros(npft,1);
lai_vp=zeros(npft,length(zzbot));

isi_a = pysi_id(ipy);
isi_z = isi_a+pysi_n(ipy)-1;

npatch = sum(sipa_n(isi_a:isi_z));
lai_pa = zeros(npatch,1);



ip=0;
for isi=isi_a:isi_z
	
  ipa_a = sipa_id(isi);
  ipa_z = ipa_a+sipa_n(isi)-1;
  
  for ipa=ipa_a:ipa_z
    
    ico_a = paco_id(ipa);
    ico_z = ico_a+paco_n(ipa)-1;
    
    afrac = area_pa(ipa)*area_si(isi);

    if(paco_n(ipa)>0)
    agb_pa(ipa) = sum(agb_co(ico_a:ico_z).*nplant(ico_a:ico_z));

ip=ip+1;
lai_pa(ip) = sum(lai_co(ico_a:ico_z));
    
    for ico=ico_a:ico_z
      
      % Find the location in the profile where this cohorts
      % leaves start and stop
      ipft = pft_co(ico);
      
      lai_pft(ipft)=lai_pft(ipft)+lai_co(ico)*afrac;
      
      hbotcrown = h2crownbh(height(ico),pft_co(ico));
      htopcrown = height(ico);
      
      [jnk kapv] = sort((hbotcrown-zzbot(:))...
			.*(zztop(:)-hbotcrown),1,'descend');
      kap=kapv(1);
      if(abs(jnk(1))<0.00000001)
	kaf = kap;
      else
	kaf = kap+1;
      end
      
      [jnk kzfv] = sort((htopcrown-zzbot(:))...
			.*(zztop(:)-htopcrown),1,'descend');
      kzp=kzfv(1);
      if(abs(jnk(1))<0.000000001)
	kzf=kzp;
      else
	kzf       = kzp-1;
      end
      
      f_depth = htopcrown-hbotcrown;  % full depth of crown
      s_depth = 0.0;                  % depth sum counter
      
      % Add some leaf area to the first partial layer
      if kap~=kaf && kap~=kzp
	p_depth     = zztop(kap)-hbotcrown;
	if(p_depth<0);display(p_depth);pause;end;
	s_depth     = s_depth+p_depth;
	lai_vp(ipft,kap) = lai_vp(ipft,kap)+lai_co(ico)*(p_depth/f_depth)*afrac;
      end
      
      % Add some leaf area to all full layers
      if(kzf>=kaf)
	for k=kaf:kzf
	  p_depth     = min([zztop(k)-zzbot(k) f_depth]);
	  if(p_depth<0);display(p_depth);pause;end;
	  s_depth     = s_depth+p_depth;
	  lai_vp(ipft,k) = lai_vp(ipft,k)+lai_co(ico)*(p_depth/f_depth)*afrac;
	end
      end
      
      % Add some leaf area to the last partial layer
      if kzp~=kzf && kzp~=kap
	p_depth     = htopcrown-zzbot(kzp);
	if(p_depth<0);display(p_depth);pause;end;
	s_depth     = s_depth+p_depth;
	lai_vp(ipft,kzp) = lai_vp(ipft,kzp)+lai_co(ico)*(p_depth/f_depth)*afrac;
      end
      
      % Final case for where the partial layers are the same
      if kzp==kap
	p_depth = f_depth;
	s_depth = s_depth+p_depth;
	lai_vp(ipft,kzp) = lai_vp(ipft,kzp)+lai_co(ico)*(p_depth/f_depth)*afrac;
      end
      
      if(abs(s_depth-f_depth)>0.001)
	display(sprintf(...
	    'DEPTH ISSUE: s_depth %d f_depth %d\n',s_depth, ...
	    f_depth));
	pause;
      end
      
    end   % for ico
      end % if(paco_n>0)
    
  end
end



% Create a figure portrait figure (for document)
% Allow for 3/4" margins = 7x9.5

figh=figure('visible',visible);
set(gcf,'PaperPositionMode','manual','Units','inches');
set(gcf,'Position',[1 1 6.0 6.0]);


rectangle('Position',[0.0 0.0 1.0 1.0],'FaceColor',[0.7 0.7 0.7]);
subplot('Position',[0.07 0.9 0.25 0.1]);
axis off;
text(0.0,0.6,...
     sprintf('%s:  %5.2fN %5.2fE',poistr,lat(ipy),lon(ipy)),...
     'FontSize',fasz+1);

npatch = length(area_pa);

for ipa=1:npatch
  pstr{ipa} = '';  %WE DONT NEED LABELS
end

comap = flipud(colormap('Summer'));
colormap(comap);

% these pie plots look funny when multiple sites are used
% and you dont sort

[age_srt,id_srt]=sort(age_pa,'descend');
area_srt = area_pa(id_srt);
agb_srt = agb_pa(id_srt);
lai_srt = lai_pa(id_srt);
dist_srt = dist_pa(id_srt);


%  PATCH AGE  (UL)
subplot('Position',[0.07 0.63 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch Age [years]','FontSize',fasz-1);
setpiecolor(hp,npatch,age_srt,comap,fasz-1);
colorbar('SouthOutside');


%  PATCH LAI  (UM)
subplot('Position',[0.41 0.67 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch LAI [m2/m2]');
setpiecolor(hp,npatch,lai_srt,comap,fasz-1);
colorbar('SouthOutside');

% PATCH AGB   (UR)
subplot('Position',[0.7 0.67 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch AGB [kg/m2]');
setpiecolor(hp,npatch,agb_srt,comap,fasz-1);
colorbar('SouthOutside');

% PATCH DISTURBANCE
subplot('Position',[0.50 0.40 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(agb_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Disturbance Regime ');
nco=length(comap);
minc = 1;
maxc = 6;
caxis([minc maxc]);
for ipa=1:npatch
id = round((nco-1)*(dist_srt(ipa)-minc)/(maxc-minc))+1;
comat(ipa,:) = comap(id,:);
end
for ip=1:npatch
 ip1 = ip*2-1;
 ip2 = ip*2;
 set(hp(ip1),'FaceColor',comat(ip,:));
 set(hp(ip2),'FontSize',fasz-1);
end
colorbar('YTick',[1 2 3 4 5 6],'YTickLabel',...
	 {'Pasture','Plantation','Treefall','Burnt','Abandoned','Logged'},...
	 'FontSize',fasz-1);%,'Location','East');


% LAI Profile (LR)
subplot('Position',[0.1 0.12 0.35 0.35]);
set(gca,'FontSize',fasz-1);
use_pft=find(lai_pft)';

hold on;
for k=1:length(zzbot)
  xloc=0;
  lai_vp(:,k) = lai_vp(:,k)./(zztop(k)-zzbot(k));						  
  for ipft=use_pft
    h1=rectangle('Position',[xloc,zzbot(k),lai_vp(ipft,k)+1e-5, ...
		    zztop(k)-zzbot(k)],'FaceColor',greeniness(ipft,:),'EdgeColor', ...
		 'k');
    xloc=xloc+lai_vp(ipft,k)+1e-5;
  end
end
xlabel('LAI [m2/m3]','FontSize',fasz-1);
ylabel('Elevation [m]','FontSize',fasz-1);
ylim([0 ymax]);
xlim([0 max([0.25 1.1*max(sum(lai_vp,1)) ])]);
grid on;
box on;
hold off;


%AGB DISTRUBUTION

bdbh = linspace(10,ndbh*10,ndbh);
for i=1:ndbh-1
  sdbh{i} = num2str(bdbh(i));
end
sdbh{ndbh}='+';
p_ids  = find(use_pft<5 & use_pft>1);
p_pfts = use_pft(p_ids);
for i=1:length(p_ids)
  spft{i} = sprintf('%s',pftshort{p_pfts(i)});
end

agb_rawt = zeros(length(p_ids),ndbh);
for i=1:length(p_pfts)
     agb_rawt(i,:) = agb_raw(use_pft(p_ids(i)),:);
end

if(prod(size(agb_rawt))>0)

subplot('Position',[0.54 0.12 0.37 0.25]);
set(gca,'FontSize',fasz-1);
bar(agb_rawt','stack');
set(gca,'xticklabel',sdbh,'FontSize',fasz-1);
grid on;
box on;
xlabel('DBH [cm]','FontSize',fasz-1);
xlim([0 12])
%view(-65,55);
title('AGB [kg/m2]','FontSize',fasz-1);
legend(spft,'Location','NorthWest');
end


% LEGEND FOR THE LAI PROFILE
ax3 = axes('Position',[0.1 0.46 0.2 0.10]);
xlim([0 1]);
ylim([0 1]);
axis off;
yloc=1.0;
for ipft=use_pft
  yloc=yloc-0.25;
  rectangle('Position',[0.05 yloc 0.1 0.25],...
	    'FaceColor',greeniness(ipft,:));
  
  tlstr=sprintf('%s (%4.2f)', pftname{ipft},lai_pft(ipft));
  text(0.20,yloc+0.05,tlstr,'FontSize',fasz-1);
  
end


oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-dpng', sprintf('%s.png',figname), '-r300');
print('-depsc', sprintf('%s.eps',figname), '-r300');
drawnow;
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)



%set(gcf,'PaperPositionMode','auto');
%print('-depsc',sprintf('%s_profbar.eps',poistr));


%print('-dpng','-r300',figname);

