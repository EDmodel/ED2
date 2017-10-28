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
           'Early Savannah',     ...   %12
           'Mid Savannah',        ...   %13
           'Late Savannah',     ...   %14
           'Araucaria',        ...   %15
           'C3 Tropical',    ...   %16
           'Liana'};      %17

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
           'Early Sav.',     ...   %12
           'Mid Sav.',        ...   %13
           'Late Sav.',     ...   %14
           'Araucaria',        ...   %15
           'C3 Trop.',    ...   %16
           'Liana'};      %17


pftshort = {'C4G',       ...   %1
           'ETR', ...   %2
           'MTR',   ...   %3
           'LTR',  ...   %4
           'TTG',   ...   %5
           'NPN',  ...   %6
           'SPN',  ...   %7
           'LCN',   ...   %8
           'EHW', ...   %9
           'MHW.',   ...   %10
           'LHW',  ...   %11
           'ESV',     ...   %12
           'MSV',        ...   %13
           'LSV',     ...   %14
           'ARC',        ...   %15
           'C3G',    ...   %16
           'LIA'};      %17



pftcolor = [221   204  119; ...
            131   204  192; ...
             68   170  153; ...
             24   102   89; ...
            235   224  170; ...
            136   204  238; ...
            182   224  245; ...
             49   112  143; ...
            204   131  192; ...
            170    68  153; ...
            102    24   89; ...
            204   131  155; ...
            170    68  102; ...
            102    24   50; ...
            115   101  184; ...
            133   118   43; ...
             51    34  136]./255;
