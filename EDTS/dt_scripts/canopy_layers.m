function [zzbot zztop] = canopy_layers(zmin,zmax,hmax,stretch);


%zmin: width of smallest (first) layer
%zmax: maximum allowable layer width
%hmax: what is you highest tree?
%stretch: what kind of stretching per layer do you want?

ztop = 0;
zbot = 0;
ilyr = 0;
while(zbot<hmax)
  ilyr   = ilyr+1;
  dz     = min([zmin*stretch.^(ilyr-1) zmax]);
  zbot   = ztop;
  ztop   = ztop+dz;
end

ncanlyr   = ilyr-1;

zztop=zeros(ncanlyr,1);
zzbot=zeros(ncanlyr,1);

zztop(1) = zmin;
zzbot(1) = 0.0;

for ilyr =2:ncanlyr
  dz           = min([zmin*stretch.^(ilyr-1) zmax]);
  zztop (ilyr) = zztop(ilyr-1)+dz;
  zzbot (ilyr) = zztop(ilyr-1);
end

end