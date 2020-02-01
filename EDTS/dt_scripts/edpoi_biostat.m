function edpoi_biostat(filename,ipy,poistr,figname,status,visible,type)

global fasz;
global pftcolor;
global pftshort;
global dtycolor;
global dtyshort;
global npft;
global npftmx;
global ndty;


% ========================================================
%
% Provide an ed state file. Define which polygon to use.
% Also, if this is part of a set, give it a name poistr.


% Get the vertical canopy layering
zmin = 0.40;
zmax = 1.0;
hmax = 42.0;
stretch=1.1;
[zzbot zztop] = canopy_layers(zmin,zmax,hmax,stretch);

ymax = hmax+2;
xmax = 0.5;



ddbh = 10;
ndbh = 11;
bdbh = linspace(ddbh,ndbh*ddbh,ndbh);

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
dbh         = hdf5read(filename,'/DBH');
pft_co      = hdf5read(filename,'/PFT');
lai_co      = hdf5read(filename,'/LAI_CO');
nplant      = hdf5read(filename,'/NPLANT');
end



agb_pa = zeros(sum(sipa_n),1);
area_dty = zeros(ndty,1);

% Process data

lai_pft = zeros(npft,1);
lai_vp=zeros(npft,length(zzbot));

isi_a = pysi_id(ipy);
isi_z = isi_a+pysi_n(ipy)-1;

npatch = sum(sipa_n(isi_a:isi_z));
lai_pa = zeros(npatch,1);

agb_sz=zeros(npft,ndbh);
agb_pft=zeros(npft,1);


ipg=0;
for isi=isi_a:isi_z
  ipa_a = sipa_id(isi);
  ipa_z = ipa_a+sipa_n(isi)-1;
  
  for ipa=ipa_a:ipa_z
    ipg=ipg+1;
    
    ico_a = paco_id(ipa);
    ico_z = ico_a+paco_n(ipa)-1;
    
    afrac = area_pa(ipa)*area_si(isi);
    idt = dist_pa(ipa);
    area_dty(idt) = area_dty(idt) + afrac;

    if(paco_n(ipa)>0)
    agb_pa(ipg) = sum(agb_co(ico_a:ico_z).*nplant(ico_a:ico_z));

    lai_pa(ipg) = sum(lai_co(ico_a:ico_z));
    
    for ico=ico_a:ico_z
      
      % Find the location in the profile where this cohorts
      % leaves start and stop
      ipft = pft_co(ico);
      lai_pft(ipft)=lai_pft(ipft)+lai_co(ico)*afrac;
      hbotcrown = h2crownbh(height(ico),pft_co(ico));
      htopcrown = height(ico);

      % Integrate size.
      idbh = max(1,min(ceil(dbh(ico)./ddbh),ndbh));
      agb_sz(ipft,idbh) = agb_sz(ipft,idbh) + nplant(ico)*agb_co(ico)*afrac;
      agb_pft(ipft)     = agb_pft(ipft)     + nplant(ico)*agb_co(ico)*afrac;


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

      if (f_depth > 0)

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
     sprintf('%s (%s):  %5.2fN %5.2fE',poistr,status,lat(ipy),lon(ipy)),...
     'FontSize',fasz+1);

npatch = length(area_pa);

for ipa=1:npatch
  pstr{ipa} = '';  %WE DONT NEED LABELS
end

comap = [140    81   10; ...
         191   129   45; ...
         223   194  125; ...
         246   232  195; ...
         235   245  245; ...
         199   234  229; ...
         128   205  193; ...
          53   151  143; ...
           1   102   94]./255;

colormap(comap);

% these pie plots look funny when multiple sites are used
% and you dont sort

[age_srt,id_srt]=sort(age_pa,'descend');
area_srt = area_pa(id_srt);
agb_srt = agb_pa(id_srt);
lai_srt = lai_pa(id_srt);
dist_srt = dist_pa(id_srt);


%  PATCH AGE  (UL)
subplot('Position',[0.05 0.65 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch Age [years]','FontSize',fasz-1);
setpiecolor(hp,npatch,age_srt,comap,fasz-1);
colorbar('EastOutside');


%  PATCH LAI  (UM)
subplot('Position',[0.35 0.65 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch LAI [m2/m2]');
setpiecolor(hp,npatch,lai_srt,comap,fasz-1);
colorbar('EastOutside');

% PATCH AGB   (UR)
subplot('Position',[0.65 0.65 0.25 0.25]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(area_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Patch AGB [kg/m2]');
setpiecolor(hp,npatch,agb_srt,comap,fasz-1);
colorbar('EastOutside');

% PATCH DISTURBANCE
subplot('Position',[0.35 0.45 0.15 0.15]);
set(gca,'FontSize',fasz-1);
paxp = zeros(size(agb_pa));
hp=pie(double(area_srt),paxp,pstr);
title('Disturbance Regime ');
for ipa=1:npatch
   ilu = dist_srt(ipa);
   comat(ipa,:) = dtycolor(ilu,:);
end
for ip=1:npatch
 ip1 = ip*2-1;
 ip2 = ip*2;
 set(hp(ip1),'FaceColor',comat(ip,:));
 set(hp(ip2),'FontSize',fasz-1);
end
colorbar('off');

% LAI Profile (LR)
subplot('Position',[0.10 0.10 0.35 0.30]);
set(gca,'FontSize',fasz-1);
use_pft=find(lai_pft)';

hold on;
for k=1:length(zzbot)
  xleft=0;
  ybottom=zzbot(k);
  delta_y=zztop(k)-zzbot(k);
  lai_vp(:,k) = lai_vp(:,k)./(zztop(k)-zzbot(k));						  
  for ipft=use_pft
    delta_x=lai_vp(ipft,k)+1e-5;
    h1=rectangle('Position' ,[xleft,ybottom,delta_x,delta_y], ...
		 'FaceColor',pftcolor(ipft,:),...
                 'EdgeColor',pftcolor(ipft,:));
    xleft=xleft+delta_x;
  end
end
xlabel('LAD [m2/m3]','FontSize',fasz-1);
ylabel('Height [m]','FontSize',fasz-1);
ylim([0 ymax]);
xlim([0 max([0.25 1.1*max(sum(lai_vp,1)) ])]);
grid on;
box on;
hold off;


%AGB DISTRUBUTION
subplot('Position',[0.55 0.10 0.35 0.30]);
set(gca,'FontSize',fasz-1);
use_pft=find(agb_pft)';
hold on;
for idbh=1:ndbh
   ybottom=0;
   xleft=bdbh(idbh)-ddbh;
   for ipft=use_pft
      delta_y=agb_sz(ipft,idbh)+1e-5;
      v1=rectangle('Position' ,[xleft,ybottom,ddbh,delta_y],...
                   'FaceColor',pftcolor(ipft,:),...
                   'EdgeColor',pftcolor(ipft,:));
      ybottom=ybottom+delta_y;
   end
end
xlabel('DBH [cm]','FontSize',fasz-1);
ylabel('AGB [kgC/m2]','FontSize',fasz-1);
xlim([0 bdbh(ndbh)]);
ylim([0 max([0.5 1.1*max(sum(agb_sz,1)) ])]);
grid on;
box on;
hold off;



% LEGEND FOR THE LAI PROFILE
ax3 = axes('Position',[0.10 0.45 0.20 0.20]);
xlim([0 1]);
ylim([0 1]);
axis off;
dy=1.0 ./ npftmx;
yloc=0;
for ipft=use_pft
  rectangle('Position',[0.05 yloc 0.1 dy],'FaceColor',pftcolor(ipft,:));

  tlstr=sprintf('%s (%4.2f)', pftshort{ipft},lai_pft(ipft));
  text(0.20,yloc + 0.5 .* dy,tlstr,'FontSize',fasz-3,...
       'HorizontalAlignment','left','VerticalAlignment','middle');
  yloc=yloc+dy;
end



% LEGEND FOR THE AGB PROFILE
ax3 = axes('Position',[0.75 0.45 0.20 0.20]);
axis off;
xlim([0 1]);
ylim([0 1]);
dy=1.0 ./ npftmx;
yloc=0;
for ipft=use_pft
  rectangle('Position',[0.05 yloc 0.1 dy],...
            'FaceColor',pftcolor(ipft,:));
  
  tlstr=sprintf('%s (%4.2f)', pftshort{ipft},agb_pft(ipft));
  text(0.20,yloc + 0.5 .* dy,tlstr,'FontSize',fasz-3,...
       'HorizontalAlignment','left','VerticalAlignment','middle');
  yloc=yloc+dy;
end


% LEGEND FOR THE DTY
use_dty = find(area_dty)';
ax3 = axes('Position',[0.55 0.45 0.20 0.175]);
xlim([0 1]);
ylim([0 1]);
axis off;
yoff=1.0 ./ ndty;
yloc=0;
for idty=use_dty
  rectangle('Position',[0.05 yloc 0.10 yoff],'FaceColor',dtycolor(idty,:));
  text(0.20,yloc + 0.5 .* dy,dtyshort(idty),'FontSize',fasz-1,...
       'HorizontalAlignment','left','VerticalAlignment','middle');
  yloc=yloc+yoff;
end


oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches','PaperPosition',newpos)
print('-dpng', sprintf('%s.png',figname), '-r300');
print('-depsc', sprintf('%s.eps',figname), '-r300');
drawnow;
set(gcf,'Units',oldscreenunits,'PaperUnits',oldpaperunits,...
    'PaperPosition',oldpaperpos)



%set(gcf,'PaperPositionMode','auto');
%print('-depsc',sprintf('%s_profbar.eps',poistr));


%print('-dpng','-r300',figname);

