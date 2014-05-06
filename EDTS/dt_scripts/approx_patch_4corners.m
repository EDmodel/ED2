function [lonpcrns,latpcrns] = ...
    approx_patch_4corners(latitude,longitude)

  % =======================================================================
  % Create the bounding boxes for each of the cells
  % This script assumes a completely amorphous geometry.
  % =======================================================================
  

  npolygons = length(latitude);
  
  lonpcrns = zeros(4,npolygons);
  latpcrns = zeros(4,npolygons);

  nlonsolved = 0;dlon_ave=0;
  nlatsolved = 0;dlat_ave=0;

  dlon_sav = zeros(npolygons,1);
  dlat_sav = zeros(npolygons,1);

  rad2s=zeros(npolygons,1);

  % Get the average radius
  for ipy = 1:npolygons
    
    radii = sqrt((latitude(ipy)-latitude).^2+(longitude(ipy)- ...
					      longitude).^2);
    [orad,irad] = sort(radii);
    
    rad2s(ipy)     = orad(2);
  end  

  rad2=mean(rad2s);
  smallrad = 0.85*(sqrt(2.0)/2.0)*rad2;
  
  for ipy = 1:npolygons
      
      %Find closest

%      display(sprintf('bounding polygon: %d of %d\n',ipy, ...
%		      npolygons));
%      display(sprintf('@ lat: %f  lon: %f\n',latitude(ipy),longitude(ipy)));
      
      
      radii = sqrt((latitude(ipy)-latitude).^2+(longitude(ipy)- ...
						longitude).^2);
      [orad,irad] = sort(radii);
      
      
      % Find the difference vectors
      
      dlon_vec = abs(longitude(ipy)-longitude);
      dlat_vec = abs(latitude(ipy)-latitude);

      % Reduce the difference vectors to those polgons within 2  radius
      ids         = find(dlon_vec<2.*rad2);
      
      if(length(ids)<1) 
	display('fail-0');
	keyboard;
      end
      
      olon_temp   = dlon_vec(ids);
      olat_temp   = dlat_vec(ids);
      
      ids         = find(olat_temp<2.*rad2);
      
      if(length(ids)<1)
	display('fail-1');
	keyboard;
	pause;
      end
      
      flats       = olat_temp(ids);
      flons       = olon_temp(ids);

      
      % Sort the latitude vector so that it descends, find the
      % closest longitude that is not on top of itself
      % Since we are going in descending order, the last is the best
            
      [slats,ilats] = sort(flats','descend');
      
      dlon=180.;
      for i=ilats
	if(flons(i)>smallrad && flons(i)<dlon)
	  dlon = flons(i);
	end
      end

      % If there are no elegible delta-lons, then we sort
      % the polygons in terms of radial distance, and then
      % recursivly find and answer

      if(dlon>179)
	dlon_solved(ipy)=0;
	display(sprintf('dlon 2 %4.3f\n',dlon));
      else
%	display(sprintf('dlon 1 %4.3f\n',dlon));
	dlon_solved(ipy)=1;
	dlon_ave=dlon_ave+dlon;
	nlonsolved=nlonsolved+1;
	dlon_sav(ipy)=dlon;	
      end
      
      [slons,ilons] = sort(flons','descend');
      dlat=180.;
      for i=ilons
	if(flats(i)>smallrad && flats(i)<dlat)
	  dlat = flats(i);
	end
      end
      
      if(dlat>179)
	dlat_solved(ipy)=0;
      else
	dlat_solved(ipy)=1;
	dlat_ave=dlat_ave+dlat;
	nlatsolved=nlatsolved+1;
	dlat_sav(ipy)=dlat;
      end

      
      %UL
      lonpcrns(1,ipy) = longitude(ipy) - 0.55*dlon;
      latpcrns(1,ipy) = latitude(ipy) + 0.55*dlat;
      
      %UR
      lonpcrns(2,ipy) = longitude(ipy) + 0.55*dlon;
      latpcrns(2,ipy) = latitude(ipy) + 0.55*dlat;
      
      %LR
      lonpcrns(3,ipy) = longitude(ipy) + 0.55*dlon;
      latpcrns(3,ipy) = latitude(ipy) - 0.55*dlat;
      
      %LL
      lonpcrns(4,ipy) = longitude(ipy) - 0.55*dlon;
      latpcrns(4,ipy) = latitude(ipy) - 0.55*dlat;
      
end      


dlat_ave=dlat_ave./nlatsolved;
dlon_ave=dlon_ave./nlonsolved;

% Check to see which ones were not solved and copy closest neighbor
% data

[ids]=find(dlon_solved==0);

for i=1:length(ids)
  ipy=ids(i);
  
  radii = sqrt((latitude(ipy)-latitude).^2+(longitude(ipy)- ...
					    longitude).^2);
  [orad,irad] = sort(radii);

  for j=2:length(irad)
    jpy=irad(j);
    if(dlon_solved(jpy)==1)
      dlon=abs(lonpcrns(2,jpy)-longitude(jpy));
      break;
    end
  end
  dlon_sav(ipy)=dlon_ave;
  dlon = dlon_ave;

  lonpcrns(1,ipy) = longitude(ipy) - 0.55*dlon;
  lonpcrns(2,ipy) = longitude(ipy) + 0.55*dlon;
  lonpcrns(3,ipy) = longitude(ipy) + 0.55*dlon;
  lonpcrns(4,ipy) = longitude(ipy) - 0.55*dlon;
end

[ids]=find(dlat_solved==0);

for i=1:length(ids)
  ipy=ids(i);

  radii = sqrt((latitude(ipy)-latitude).^2+(longitude(ipy)- ...
                                            longitude).^2);
  [orad,irad] = sort(radii);

  for j=2:length(irad)
    jpy=irad(j);
    if(dlat_solved(jpy)==1)
      dlat=abs(latpcrns(1,jpy)-latitude(jpy));
      break;
    end
  end
  dlat_sav(ipy)=dlat_ave;
  dlat=dlat_ave;

  latpcrns(1,ipy) = latitude(ipy) + 0.55*dlat;
  latpcrns(2,ipy) = latitude(ipy) + 0.55*dlat;
  latpcrns(3,ipy) = latitude(ipy) - 0.55*dlat;
  latpcrns(4,ipy) = latitude(ipy) - 0.55*dlat;

end

display('COMPLETED BOUNDING RECTANGULAR POLYGONS');
