function [] = setpiecolor(hp,np,cvar,comap,fasz)

nco=length(comap);
minc = min(cvar);
maxc = max(cvar);

if(maxc>minc)
  caxis([minc maxc]);
 else
   minc=minc-0.25*(abs(minc));
   maxc=maxc+0.25*(abs(maxc));
   if(minc==0 || maxc==0)
     minc = -1;
     maxc = 1;
   end
   caxis([minc maxc]);
end

for ip=1:np
  id = round((nco-1)*(cvar(ip)-minc)/(maxc-minc))+1;
comat(ip,:) = comap(id,:);
end

for ip=1:np

 ip1 = ip*2-1;
 ip2 = ip*2;
 set(hp(ip1),'FaceColor',comat(ip,:));
 set(hp(ip2),'FontSize',fasz);
end
