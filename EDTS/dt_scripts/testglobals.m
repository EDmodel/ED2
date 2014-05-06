% Global variables


global pftcolor;
global pftshort;
global pftlong;
global pftname;

pftlong = {'C4 Grass',       ...   %1
           'Early Tropical', ...   %2
           'Mid Tropical',   ...   %3
           'Late Tropical',  ...   %4
           'C3 Temperate',   ...   %5
           'North Pine',  ...   %6
           'South Pine',  ...   %7
           'Late Conifer',   ...   %8
           'Early Hardwood', ...   %9
           'Mid Hardwood',   ...   %10
           'Late Hardwood',  ...   %11
           'C3 Pasture',     ...   %12
           'C3 Crop',        ...   %13
           'C4 Pasture',     ...   %14
           'C4 Crop',        ...   %15
           'C3 Tropical',    ...   %16
           'Araucaria Pine'};      %17

pftname = {'C4 Grass',       ...   %1
           'Early Trop.', ...   %2
           'Mid Trop.',   ...   %3
           'Late Trop.',  ...   %4
           'C3 Temp.',   ...   %5
           'North Pine',  ...   %6
           'South Pine',  ...   %7
           'Late Conifer',   ...   %8
           'Early Hard.', ...   %9
           'Mid Hard.',   ...   %10
           'Late Hard.',  ...   %11
           'C3 Past.',     ...   %12
           'C3 Crop',        ...   %13
           'C4 Past.',     ...   %14
           'C4 Crop',        ...   %15
           'C3 Trop.',    ...   %16
           'Araucaria Pine'};      %17


pftshort = {'C4',       ...   %1
           'ETROP', ...   %2
           'MTROP',   ...   %3
           'LTROP',  ...   %4
           'C3TEMP',   ...   %5
           'NPINE',  ...   %6
           'SPINE',  ...   %7
           'LCON',   ...   %8
           'EHARD', ...   %9
           'MHARD.',   ...   %10
           'LHARD',  ...   %11
           'C3PAST',     ...   %12
           'C3CROP',        ...   %13
           'C4PAST',     ...   %14
           'C4CROP',        ...   %15
           'C3TROP',    ...   %16
           'APINE'};      %17



pftcolor = [255   215    0; ...
            127   255    0; ...
             69   139    0; ...
              0    78    0; ...
            171   130  255; ...
              0   191  255; ...
             72   209  204; ...
             39    64  139; ...
	    255   140    0; ...
	    255    69    0; ...
	    178    34   34; ...
	     85    26  139; ...
	    191    62  255; ...
	    184   134   11; ...
	    240   230  140; ...
	    205   190  112; ...
             79   148  205]./255;
