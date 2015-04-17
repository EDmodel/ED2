%==========================================================================
% Program: dev_test.m
%
% This script runs a suit of checks on EDM output.  The purpose
% of the checks is to compare newly developed model runs to the
% the mainline.
%
% Please read the README file for more information
%
%==========================================================================

clear all;
close all;

%==========================================================================
%     User defined variables
%==========================================================================

test_name = 'r85ghubrapid';

use_m34 = true;       % SOI Manaus km34
use_ata = true;       % SOI Atacama
use_s67 = true;       % SOI Santarem km 67
use_har = true;       % SOI Harvard Forest
use_pdg = true;       % SOI Pe de Gigante
use_cax = true;       % SOI Caxuana
use_ton = true;       % SOI Tonzi (temperate)
use_tnf = true;       % SOI Tapajos National Forest
use_pet = true;       % SOI Petrolina
use_hip = true;       % SOI Petrolina (short high frequency)
use_him = true;       % SOI Manaus (short high frequency)
use_rjg = true;       % GRIDDED centered on Rio Jaru

%==========================================================================



site_name  = {'Manaus km 34', ...
              'Atacama Desert',...
              'Santarem km 67',...
              'Harvard Forest',...
              'Pe de Gigante',...
              'Caxuana',...
              'Tonzi',...
              'Tapajos National Forest',...
              'Petrolina'};

siteid     = {'m34',...
              'ata',...
              's67',...
              'har',...
              'pdg',...
              'cax',...
              'ton',...
              'tnf',...
              'pet'};

hifr_name = {'Petrolina High Frequency','Manaus High Frequency'};
hifrid    = {'hip','him'};


gridid     = {'rjg'};
grid_name  = {'12x12 Offline Grid - Rebio Jaru'};

addpath(strcat(pwd,'/dt_scripts'));
addpath(strcat(pwd,'/dt_scripts/cbfreeze'));
addpath(strcat(pwd,'/dt_scripts/exportfig'));

%==========================================================================

%set(0,'DefaultAxesFontName','Courier 10 Pitch');
global fasz;
fasz = 10;
visible = 'off';
addpath('dt_scripts');
testglobals;
outdir = strcat(test_name,'/report/');
mkdir(outdir);


nsite = numel(siteid);
nhifr = numel(hifrid);
ngrid = numel(gridid);

%==========================================================================
% Read in the xml data to help write the report
%==========================================================================

xmlfile = strcat(test_name,'/test_text.xml');
try
    xmltree = xmlread(xmlfile);
catch
    error('Failed to read XML file: %s',xmlfile);
end

xmlitems = xmltree.getElementsByTagName('description');
xmlitem0 = xmlitems.item(0);
branch_item = xmlitem0.getElementsByTagName('branch_version');
branch_version = char(branch_item.item(0).getFirstChild.getData);
committer_item = xmlitem0.getElementsByTagName('committer_name');
committer_name = char(committer_item.item(0).getFirstChild.getData);
tester_item = xmlitem0.getElementsByTagName('tester_name');
tester_name = char(tester_item.item(0).getFirstChild.getData);
test_item          = xmlitem0.getElementsByTagName('test_description');
test_description   = char(test_item.item(0).getFirstChild.getData);


%==========================================================================
%  Determine the status of the runs
%  If not *_out files exist in the report, then don't use this one
%  Use the check script? COMP-FAILED-RUNNING-DNEXIST
%==========================================================================

display(sprintf('\nThe following sites will be assessed:\n'));

use_site = [use_m34,use_ata,use_s67,...
            use_har,use_pdg,use_cax,...
            use_ton,use_tnf,use_pet];
        
use_hifr = [use_hip,use_him];

use_grid = [use_rjg];

for is=1:nsite
    testout_srch = sprintf('%s/test_%s.',test_name, ...
        siteid{is});
    mainout_srch  = sprintf('%s/main_%s.',test_name, ...
        siteid{is});
    dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
        siteid{is});

    srch=dir(strcat(testout_srch,'*out'));
    testout_str=sprintf('%s/%s',test_name,srch(end).name);

    srch=dir(strcat(dbugout_srch,'*out'));
    dbugout_str=sprintf('%s/%s',test_name,srch(end).name);

    srch=dir(strcat(mainout_srch,'*out'));
    mainout_str=sprintf('%s/%s',test_name,srch(end).name);

    if(use_site(is))
    if (exist(testout_str,'file') && ...
        exist(mainout_str,'file') && ...
        exist(dbugout_str,'file'))
        use_site(is)=true;
        display(sprintf('%s - %s',siteid{is},site_name{is}));
    else
        use_site(is)=false;
    end
    end
end

for ih=1:nhifr
    testout_srch = sprintf('%s/test_%s.',test_name, ...
        hifrid{ih});
    mainout_srch  = sprintf('%s/main_%s.',test_name, ...
        hifrid{ih});
    dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
        hifrid{ih});

    srch=dir(strcat(testout_srch,'*out'));
    testout_str=sprintf('%s/%s',test_name,srch(end).name);

    srch=dir(strcat(dbugout_srch,'*out'));
    dbugout_str=sprintf('%s/%s',test_name,srch(end).name);

    srch=dir(strcat(mainout_srch,'*out'));
    mainout_str=sprintf('%s/%s',test_name,srch(end).name);            

    if(use_hifr(ih))
    if (exist(testout_str,'file') && ...
        exist(mainout_str,'file') && ...
        exist(dbugout_str,'file'))
        use_hifr(ih)=true;
        display(sprintf('%s - %s',hifrid{ih},hifr_name{ih}));
    else
        use_hifr(ih)=false;
    end
    end
end

for ig=1:ngrid

    testout_srch = sprintf('%s/test_%s.',test_name, ...
        gridid{ig});
    mainout_srch  = sprintf('%s/main_%s.',test_name, ...
        gridid{ig});
    dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
        gridid{ig});

    srch=dir(strcat(testout_srch,'*out'));
    testout_str=sprintf('%s/%s',test_name,srch(end).name);
    srch=dir(strcat(dbugout_srch,'*out'));
    dbugout_str=sprintf('%s/%s',test_name,srch(end).name);
    srch=dir(strcat(mainout_srch,'*out'));
    mainout_str=sprintf('%s/%s',test_name,srch(end).name);

    if(use_grid(ig))
    if (exist(testout_str,'file') && ...
        exist(mainout_str,'file') && ...
        exist(dbugout_str,'file'))
        use_grid(ig)=true;
        display(sprintf('%s - %s',gridid{ig},grid_name{ig}));
    else
        use_grid(ig)=false;
    end
    end
end

display(sprintf('\nThe Following sites are not available:\n'));
ii=0;
for is=1:nsite
    if(~use_site(is))
        display(sprintf('%s - %s\n',siteid{is},site_name{is}));
        ii=ii+1;
    end
end

for ih=1:nhifr
    if(~use_hifr(ih))
        display(sprintf('%s - %s\n',hifrid{ih},hifr_name{ih}));
        ii=ii+1;
    end
end

for ig=1:ngrid
    if(~use_grid(ig))
        display(sprintf('%s - %s\n',gridid{ig},grid_name{ig}));
        ii=ii+1;
    end
end

if(ii==0)
    display(sprintf('None'));
end


% Create a new list of the sites that some output was generated

nsiteid={};
ngridid={};
nhifrid={};
temp_name = site_name;
site_name = {};
is=0;
for ii=1:nsite
    if(use_site(ii))
        is=is+1;
        nsiteid{is} = siteid{ii};
        site_name{is} = temp_name{ii};
    end
end

temp_name = grid_name;
grid_name = {};
ig=0;
for ii=1:ngrid
    if(use_grid(ii))
        ig=ig+1;
        ngridid{ig} = gridid{ii};
        grid_name{ig} = temp_name{ii};
    end
end

temp_name = hifr_name;
hifr_nam  = {};
ih=0;
for ii=1:nhifr
    if(use_hifr(ii))
        ih=ih+1;
        nhifrid{ih} = hifrid{ii};
        hifr_name{ih} = temp_name{ii};
    end
end

hifrid = nhifrid;clear nhifrid;
gridid = ngridid;clear ngridid;
siteid = nsiteid;clear nsiteid;

nsite = is;
ngrid = ig;
nhifr = ih;

runstat = {'Fail','Pass'};
pause(2);


% =========================================================================
% Part 0 - Check if debug sites completed
% =========================================================================

display(sprintf('\nChecking Simulations for Completion\n'));

if nsite>0; spass=zeros(nsite,3);end
if ngrid>0; gpass=zeros(nsite,3);end
if nhifr>0; hpass=zeros(nsite,3);end

if nsite>0

    for is=1:nsite

        dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
            siteid{is});
        srch=dir(strcat(dbugout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                spass(is,1) = 1;
            end
        end
        fclose(fid);

        testout_srch  = sprintf('%s/test_%s.',test_name, ...
            siteid{is});
        srch=dir(strcat(testout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                spass(is,2) = 1;
            end
        end
        fclose(fid);        
        
        mainout_srch  = sprintf('%s/main_%s.',test_name, ...
            siteid{is});
        srch=dir(strcat(mainout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);  

        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                spass(is,3) = 1;
            end
        end
        fclose(fid);        
        
        display(sprintf('%s - dbug:%s test:%s main:%s ',...
            siteid{is},runstat{spass(is,1)+1}, ...
            runstat{spass(is,2)+1}, ...
            runstat{spass(is,3)+1}));
        
    end
end

if ngrid>0
    for ig=1:ngrid

        dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
            gridid{ig});
        srch=dir(strcat(dbugout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                gpass(ig,1)=1;
            end
        end
        fclose(fid);
        
        testout_srch  = sprintf('%s/test_%s.',test_name, ...
            gridid{ig});
        srch=dir(strcat(testout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                gpass(ig,2)=1;
            end
        end
        fclose(fid);
        
        mainout_srch  = sprintf('%s/main_%s.',test_name, ...
            gridid{ig});
        srch=dir(strcat(mainout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);       
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                gpass(ig,3)=1;
            end
        end
        fclose(fid);
        
        display(sprintf('%s - dbug:%s test:%s main:%s ',...
            gridid{ig},runstat{gpass(ig,1)+1}, ...
            runstat{gpass(ig,2)+1}, ...
            runstat{gpass(ig,3)+1}));
    end
end

if nhifr>0
    for ih=1:nhifr
        dbugout_srch  = sprintf('%s/dbug_%s.',test_name, ...
            hifrid{ih});
        srch=dir(strcat(dbugout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                hpass(ih,1) = 1;
            end
        end
        fclose(fid);
        
        testout_srch  = sprintf('%s/test_%s.',test_name, ...
             hifrid{ih});
        srch=dir(strcat(testout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                hpass(ih,2) = 1;
            end
        end
        fclose(fid);
        
        mainout_srch  = sprintf('%s/main_%s.',test_name, ...
             hifrid{ih});
        srch=dir(strcat(mainout_srch,'*out'));
        outfile=sprintf('%s/%s',test_name,srch(end).name);
        fid=fopen(outfile);
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if(strfind(tline,'ED-2.2 execution ends'))
                hpass(ih,3) = 1;
            end
        end
        fclose(fid);
        
        display(sprintf('%s - dbug:%s test:%s main:%s ',...
            hifrid{ih},runstat{hpass(ih,1)+1}, ...
            runstat{hpass(ih,2)+1}, ...
            runstat{hpass(ih,3)+1}));
    end
end

%==========================================================================
% PART 1: Loop through SOIs
%==========================================================================

if nsite>0

display(sprintf('\nChecking SOI Fluxes, Succession and Profiles'));

latex_ftab = zeros(9,nsite);
latex_fname = {'$\\Delta ET$','$\\Delta SHF$','$\\Delta R_{net}$','$\\Delta R_{SWU}$','$\\Delta GPP$',... 
            '$\\Delta NEP$','$\\Delta CO2_{C}$','$\\Delta \\theta_{50cm}$','$\\Delta T_L$'};
latex_funit = {'$[mm/m^2]$','$[W/m^2]$','$[W/m^2]$','$[W/m^2]$'...
            '$[kgC/m^2]$','$[kgC/m^2]$','$[ppm]$','$[m^3/m^3]$','$[^oC]$'};

for is = 1:nsite
    
    if(~(spass(is,2) && spass(is,3)))
        display(sprintf('Site: %s did not complete both main and test',...
            siteid{is}));
    else
        
        display(sprintf('%s\n',siteid{is}));
        
        test_q_pfx = sprintf('%s/F_test_%s/test_%s-Q-',test_name,siteid{is},siteid{is});
        cont_q_pfx = sprintf('%s/F_main_%s/main_%s-Q-',test_name,siteid{is},siteid{is});
        test_s_pfx = sprintf('%s/S_test_%s/test_%s-S-',test_name,siteid{is},siteid{is});
        cont_s_pfx = sprintf('%s/S_main_%s/main_%s-S-',test_name,siteid{is},siteid{is});
        
        id=strfind(test_q_pfx,'/');
        test_q_dir = test_q_pfx(1:id(end));
        
        id=strfind(test_s_pfx,'/');
        test_s_dir = test_s_pfx(1:id(end));
        
        id=strfind(cont_q_pfx,'/');
        cont_q_dir = cont_q_pfx(1:id(end));
        
        id=strfind(cont_s_pfx,'/');
        cont_s_dir = cont_s_pfx(1:id(end));
        
        test_q_flist = dir(strcat(test_q_pfx,'*h5'));
        cont_q_flist = dir(strcat(cont_q_pfx,'*h5'));
        test_s_flist = dir(strcat(test_s_pfx,'*h5'));
        cont_s_flist = dir(strcat(cont_s_pfx,'*h5'));
        
        nqfiles     = length(test_q_flist);
        nsfiles     = length(test_s_flist);
        
        if (nqfiles ~= length(cont_q_flist))
            display(sprintf('Q File lists are different lengths - %s',...
                siteid{is}));
            return;
        end
        
        if (nsfiles ~= length(cont_s_flist))
            display(sprintf('S File lists are different lengths - %s',...
                siteid{is}));
            return;
        end
        
        dnq = zeros(nqfiles,1);
        dns = zeros(nsfiles,1);
        
        
        %==================================================================
        % Comparison of Flux Variables
        %==================================================================
        
        for it=1:nqfiles
            
            tqfile = strcat(test_q_dir,test_q_flist(it).name);
            cqfile = strcat(cont_q_dir,cont_q_flist(it).name);
            
            iyear  = str2double(tqfile(end-23:end-20));
            imonth = str2double(tqfile(end-18:end-17));
            idate  = str2double(tqfile(end-15:end-14));
            ihour  = str2double(tqfile(end-12:end-11));
            iminute= str2double(tqfile(end-10:end-9));
            isecond= str2double(tqfile(end-8:end-7));
            
            dnq(it) = datenum(iyear,imonth,idate,ihour,iminute,isecond);
            
            iyear  = str2double(cqfile(end-23:end-20));
            imonth = str2double(cqfile(end-18:end-17));
            idate  = str2double(cqfile(end-15:end-14));
            ihour  = str2double(cqfile(end-12:end-11));
            iminute= str2double(cqfile(end-10:end-9));
            isecond= str2double(cqfile(end-8:end-7));
            
            if(datenum(iyear,imonth,idate,ihour,iminute,isecond)~=dnq(it))
                display('Q-File Time Mismatch, check you directories');
                display('and the prefixes');
                return;
            end
            
            if(it==nqfiles)
                iyearz=iyear;
            end
            
            if(it==1)
                iyeara=iyear;
%                nhrs = numel(hdf5read(tqfile,'/QMEAN_VAPOR_AC'));

                nhrs = numel(hdf5read(tqfile,'/QMEAN_VAPOR_AC_PY'));

                tmp = -hdf5read(tqfile,'/SLZ');
                k50 = find(tmp>0.5,1,'last');
                nz = numel(tmp);
                dz = zeros(nz,1);
                for iz=1:nz-1
                    dz(iz)=tmp(iz)-tmp(iz+1);
                end
                dz(nz) = tmp(nz); 
                slz=zeros(nz,nhrs);
                for ihr=1:nhrs
                    slz(:,ihr)=dz;
                end
                
                
                phrs = linspace(24/nhrs,24,nhrs);
                npts = nqfiles*nhrs;
                nmos    = zeros(12,1);
                
                flux_tp = zeros(npts,5);
                flux_td = zeros(nhrs,5);
                flux_tm = zeros(12,5);
                flux_cp = zeros(npts,5);
                flux_cd = zeros(nhrs,5);
                flux_cm = zeros(12,5);
                
                % FAST_SOIL_C, SLOW_SOIL_C,
                
                state_tp = zeros(npts,5);
                state_td = zeros(nhrs,5);
                state_tm = zeros(12,5);
                state_cp = zeros(npts,5);
                state_cd = zeros(nhrs,5);
                state_cm = zeros(12,5);
                
                flux_names = {'ET',...
                    'SHF',...
                    'R_{net}',...
                    'R_{SWU}',...
                    'R_{LWU}'};
                
                state_names = {'GPP',...
                    'NEP',...
                    'Canopy CO2',...
                    '50cm Soil Moisture'...
                    'Leaf Temp'};
                
                flux_units =    {'[mm/m^2/mo]',...
                    '[W/m^2]',...
                    '[W/m^2]',...
                    '[W/m^2]',...
                    '[W/m^2]'};
                
                state_units = {'[kgC/m^2/yr]',...
                    '[kgC/m^2/yr]'...
                    '[ppm]',...
                    '[m^3/m^3]',...
                    '[^oC]'};
                
            end
            
            id1 = (it-1)*nhrs+1;
            id2 = (it*nhrs);
            nmos(imonth)=nmos(imonth)+1;
            
            itmp=1;
            tmp = hdf5read(tqfile,'/QMEAN_VAPOR_AC_PY')*-86400*30;
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_VAPOR_AC_PY')*-86400*30;
            flux_cp(id1:id2,itmp) = tmp;
            flux_cd(:,itmp)      = flux_cd(:,itmp) + tmp./nqfiles;
            flux_cm(imonth,itmp) = flux_cm(imonth,itmp)+mean(tmp);
            
            itmp=2;
            tmp = -hdf5read(tqfile,'/QMEAN_SENSIBLE_AC_PY');
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            
            tmp = -hdf5read(cqfile,'/QMEAN_SENSIBLE_AC_PY');
            flux_cp(id1:id2,itmp) = tmp;
            flux_cd(:,itmp)      = flux_cd(:,itmp) + tmp./nqfiles;
            flux_cm(imonth,itmp) = flux_cm(imonth,itmp)+mean(tmp);
            
            itmp=3;
            tmp = hdf5read(tqfile,'/QMEAN_RNET_PY');
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_RNET_PY');
            flux_cp(id1:id2,itmp) = tmp;
            flux_cd(:,itmp)      = flux_cd(:,itmp) + tmp./nqfiles;
            flux_cm(imonth,itmp) = flux_cm(imonth,itmp)+mean(tmp);
            
            itmp=4;
            tmp = hdf5read(tqfile,'/QMEAN_RSHORTUP_PY');
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_RSHORTUP_PY');
            flux_cp(id1:id2,itmp) = tmp;
            flux_cd(:,itmp)      = flux_cd(:,itmp) + tmp./nqfiles;
            flux_cm(imonth,itmp) = flux_cm(imonth,itmp)+mean(tmp);
            
            itmp=5;
            tmp = hdf5read(tqfile,'/QMEAN_RLONGUP_PY');
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_RLONGUP_PY');
            flux_cp(id1:id2,itmp) = tmp;
            flux_cd(:,itmp)      = flux_cd(:,itmp) + tmp./nqfiles;
            flux_cm(imonth,itmp) = flux_cm(imonth,itmp)+mean(tmp);
            
            itmp=1;
            tmp = hdf5read(tqfile,'/QMEAN_GPP_PY');
            state_tp(id1:id2,itmp) = tmp;
            state_td(:,itmp)      = state_td(:,itmp) + tmp./nqfiles;
            state_tm(imonth,itmp) = state_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_GPP_PY');
            state_cp(id1:id2,itmp) = tmp;
            state_cd(:,itmp)      = state_cd(:,itmp) + tmp./nqfiles;
            state_cm(imonth,itmp) = state_cm(imonth,itmp)+mean(tmp);
            
            itmp=2;
            tmp = hdf5read(tqfile,'/QMEAN_NEP_PY');
            state_tp(id1:id2,itmp) = tmp;
            state_td(:,itmp)      = state_td(:,itmp) + tmp./nqfiles;
            state_tm(imonth,itmp) = state_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_NEP_PY');
            state_cp(id1:id2,itmp) = tmp;
            state_cd(:,itmp)      = state_cd(:,itmp) + tmp./nqfiles;
            state_cm(imonth,itmp) = state_cm(imonth,itmp)+mean(tmp);
            
            itmp=3;
            tmp = hdf5read(tqfile,'/QMEAN_CAN_CO2_PY');
            state_tp(id1:id2,itmp) = tmp;
            state_td(:,itmp)      = state_td(:,itmp) + tmp./nqfiles;
            state_tm(imonth,itmp) = state_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_CAN_CO2_PY');
            state_cp(id1:id2,itmp) = tmp;
            state_cd(:,itmp)      = state_cd(:,itmp) + tmp./nqfiles;
            state_cm(imonth,itmp) = state_cm(imonth,itmp)+mean(tmp);
            
            itmp=4;
            delz = sum(slz(k50:nz,1));
            tmp2d = hdf5read(tqfile,'/QMEAN_SOIL_WATER_PY');
            tmp   = sum(tmp2d(k50:nz,:).*slz(k50:nz,:),1)'./delz;
            state_tp(id1:id2,itmp) = tmp;
            state_td(:,itmp)      = state_td(:,itmp) + tmp./nqfiles;
            state_tm(imonth,itmp) = state_tm(imonth,itmp)+mean(tmp);
            
            tmp2d = hdf5read(cqfile,'/QMEAN_SOIL_WATER_PY');
            tmp   = sum(tmp2d(k50:nz,:).*slz(k50:nz,:),1)'./delz;
            state_cp(id1:id2,itmp) = tmp;
            state_cd(:,itmp)      = state_cd(:,itmp) + tmp./nqfiles;
            state_cm(imonth,itmp) = state_cm(imonth,itmp)+mean(tmp);
            
            itmp=5;
            tmp = hdf5read(tqfile,'/QMEAN_LEAF_TEMP_PY')-273.14;
            state_tp(id1:id2,itmp) = tmp;
            state_td(:,itmp)      = state_td(:,itmp) + tmp./nqfiles;
            state_tm(imonth,itmp) = state_tm(imonth,itmp)+mean(tmp);
            
            tmp = hdf5read(cqfile,'/QMEAN_LEAF_TEMP_PY')-273.14;
            state_cp(id1:id2,itmp) = tmp;
            state_cd(:,itmp)      = state_cd(:,itmp) + tmp./nqfiles;
            state_cm(imonth,itmp) = state_cm(imonth,itmp)+mean(tmp);
            
        end
        
        % Normalize means
        for imo=1:12
            flux_tm(imo,:) = flux_tm(imo,:)./nmos(imo);
            flux_cm(imo,:) = flux_cm(imo,:)./nmos(imo);
            state_tm(imo,:) = state_tm(imo,:)./nmos(imo);
            state_cm(imo,:) = state_cm(imo,:)./nmos(imo);
        end
        
        % Save means for the quick reference table
        latex_ftab(1,is) = mean(flux_tm(:,1))-mean(flux_cm(:,1));
        latex_ftab(2,is) = mean(flux_tm(:,2))-mean(flux_cm(:,2));
        latex_ftab(3,is) = mean(flux_tm(:,3))-mean(flux_cm(:,3));
        latex_ftab(4,is) = mean(flux_tm(:,4))-mean(flux_cm(:,4));
        latex_ftab(5,is) = mean(state_tm(:,1))-mean(state_cm(:,1));
        latex_ftab(6,is) = mean(state_tm(:,2))-mean(state_cm(:,2));
        latex_ftab(7,is) = mean(state_tm(:,3))-mean(state_cm(:,3));
        latex_ftab(8,is) = mean(state_tm(:,4))-mean(state_cm(:,4));
        latex_ftab(9,is) = mean(state_tm(:,5))-mean(state_cm(:,5));
        
%        latex_fname = {'$ET$','$SHF$','$R_{net}$','$R_{SWU}$','$GPP$',...
%            '$NEP$','$CO2_{C}$','$\\theta_{50cm}$','$T_L$'};
        
%        latex_funit = {'$[mm/m^2]$','$[W/m^2]$','$[W/m^2]$','$[W/m^2]$'...
%            '$[kgC/m^2]$','$[kgC/m^2]$','$[ppm]$','$[m^3/m^3]$','$[^oC]$'};
        
        
        
        fluxes_img{is} = sprintf('%sfluxes_%s.eps',outdir,siteid{is});
        states_img{is} = sprintf('%sstates_%s.eps',outdir,siteid{is});
        
        titlestr = sprintf('%s (%4i - %4i)',site_name{is},iyeara,iyearz);
        
        % Make a plot
        %==================================================================
        plot_fluxes(flux_tp,flux_cp,flux_td,flux_cd,flux_tm,flux_cm,...
            flux_names,flux_units,phrs,titlestr,fluxes_img{is},visible);
        
        plot_states(state_tp,state_cp,state_td,state_cd,state_tm,state_cm,...
            state_names,state_units,phrs,titlestr,states_img{is},visible);
        
        
        %==================================================================
        % End stage biomass SOI
        %==================================================================
        
        tsfile = strcat(test_s_dir,test_s_flist(end).name);
        csfile = strcat(cont_s_dir,cont_s_flist(end).name);
        
        testfig = sprintf('%sprofbar_test_%s',outdir,siteid{is});
        edpoi_biostat(tsfile,1,sprintf('test-%s',siteid{is}),testfig,visible,2);
        
        contfig = sprintf('%sprofbar_main_%s',outdir,siteid{is});
        edpoi_biostat(csfile,1,sprintf('main-%s',siteid{is}),contfig,visible,1);
        
        strcomp_timg{is} = sprintf('%s.eps',testfig);
        strcomp_cimg{is} = sprintf('%s.eps',contfig);

        strcomp_img{is} = sprintf('%sprofbar_%s.eps',outdir,siteid{is});
        
        % Concatinate the images
        
        [image1,map1] = imread(sprintf('%s.png',contfig));
        [image2,map2] = imread(sprintf('%s.png',testfig));
        
        iwidth = max(size(image1,2),size(image2,2));
        if size(image1,2) < iwidth
            image1(1,iwidth,1) = 0;
        end
        if size(image2,2) < iwidth
            image2(1,iwidth,1) = 0;
        end
        image3 = cat(2,image1,image2);
        bpfig = figure('Visible','off');
       
        image(image3);
%        print(strcomp_img{is},'-depsc2','-r200');
        saveas(bpfig,strcomp_img{is},'eps2c')

%        imwrite(image3,strcomp_img{is},'png');



        if(nsfiles>2)
            
            %==============================================================
            % Successional dynamics
            %==============================================================
            
            for it=1:nsfiles
                
                tsfile = strcat(test_s_dir,test_s_flist(it).name);
                csfile = strcat(cont_s_dir,cont_s_flist(it).name);
                
                iyear  = str2double(tsfile(end-23:end-20));
                imonth = str2double(tsfile(end-18:end-17));
                idate  = str2double(tsfile(end-15:end-14));
                ihour  = str2double(tsfile(end-12:end-11));
                iminute= str2double(tsfile(end-10:end-9));
                isecond= str2double(tsfile(end-8:end-7));
                
                dns(it) = datenum(iyear,imonth,idate,ihour,iminute,isecond);
                
                iyear  = str2double(csfile(end-23:end-20));
                imonth = str2double(csfile(end-18:end-17));
                idate  = str2double(csfile(end-15:end-14));
                ihour  = str2double(csfile(end-12:end-11));
                iminute= str2double(csfile(end-10:end-9));
                isecond= str2double(csfile(end-8:end-7));
                
                if(datenum(iyear,imonth,idate,ihour,iminute,isecond)~=dns(it))
                    display('S-File Time Mismatch, check you directories');
                    display('and the prefixes');
                    return;
                end
                
                if(it==nsfiles)
                    iyearz=iyear;
                end
                
                if(it==1)
                    iyeara=iyear;
                    tmp = hdf5read(tsfile,'/AGB_PY');
                    agb_t = zeros([nsfiles,size(tmp,1)]);
                    agb_c = zeros([nsfiles,size(tmp,1)]);
                    lai_t = zeros([nsfiles,size(tmp,1)]);
                    lai_c = zeros([nsfiles,size(tmp,1)]);
                end
                
                tmp = hdf5read(tsfile,'/AGB_PY');
                agb_t(it,:) = sum(tmp,2);
                
                tmp = hdf5read(csfile,'/AGB_PY');
                agb_c(it,:) = sum(tmp,2);
                
                paco_id     = hdf5read(tsfile,'/PACO_ID');
                sipa_id     = hdf5read(tsfile,'/SIPA_ID');
                pysi_id     = hdf5read(tsfile,'/PYSI_ID');
                paco_n      = hdf5read(tsfile,'/PACO_N');
                sipa_n      = hdf5read(tsfile,'/SIPA_N');
                pysi_n      = hdf5read(tsfile,'/PYSI_N');
                area_si     = hdf5read(tsfile,'/AREA_SI');
                area_pa     = hdf5read(tsfile,'/AREA');
                if(sum(paco_n)>0)
                    pft_co      = hdf5read(tsfile,'/PFT');
                    lai_co      = hdf5read(tsfile,'/LAI_CO');
                end
                
                isi_a = pysi_id(1);
                isi_z = isi_a+pysi_n(1)-1;
                for isi=isi_a:isi_z
                    ipa_a = sipa_id(isi);
                    ipa_z = ipa_a+sipa_n(isi)-1;
                    for ipa=ipa_a:ipa_z
                        
                        if(paco_n(ipa)>0)
                            ico_a = paco_id(ipa);
                            ico_z = ico_a+paco_n(ipa)-1;
                            afrac = area_pa(ipa)*area_si(isi);
                            for ico=ico_a:ico_z
                                ipft = pft_co(ico);
                                lai_t(it,ipft)=lai_t(it,ipft)+lai_co(ico)*afrac;
                            end
                        end
                    end
                end
                
                paco_id     = hdf5read(csfile,'/PACO_ID');
                sipa_id     = hdf5read(csfile,'/SIPA_ID');
                pysi_id     = hdf5read(csfile,'/PYSI_ID');
                paco_n      = hdf5read(csfile,'/PACO_N');
                sipa_n      = hdf5read(csfile,'/SIPA_N');
                pysi_n      = hdf5read(csfile,'/PYSI_N');
                area_si     = hdf5read(csfile,'/AREA_SI');
                area_pa     = hdf5read(csfile,'/AREA');
                if(sum(paco_n)>0)
                    pft_co      = hdf5read(csfile,'/PFT');
                    lai_co      = hdf5read(csfile,'/LAI_CO');
                end
                isi_a = pysi_id(1);
                isi_z = isi_a+pysi_n(1)-1;
                for isi=isi_a:isi_z
                    ipa_a = sipa_id(isi);
                    ipa_z = ipa_a+sipa_n(isi)-1;
                    for ipa=ipa_a:ipa_z
                        if(paco_n(ipa)>0)
                            ico_a = paco_id(ipa);
                            ico_z = ico_a+paco_n(ipa)-1;
                            afrac = area_pa(ipa)*area_si(isi);
                            for ico=ico_a:ico_z
                                ipft = pft_co(ico);
                                lai_c(it,ipft)=lai_c(it,ipft)+lai_co(ico)*afrac;
                            end
                        end
                    end
                end
            end
            
            pftsucc_img{is} = sprintf('%sagb_lai_pft_%s.eps',outdir,siteid{is});
            titlestr = sprintf('%s\n',site_name{is});
            
            plot_succession(dns,agb_t,agb_c,lai_t,lai_c,titlestr, ...
                pftsucc_img{is},visible)
            
            pftsucc_plt(is)=1;
        else
            pftsucc_plt(is)=0;
        end % if nsfiles>2
        
    end % if passed
    
end     % for is=1:npoi
end     % if nsite>0


if nhifr>0

display(sprintf('\nHigh Frequency Output'))
display(sprintf('Checking: Patch Level Mass and Energy Conservation\n'));    

latex_htab = zeros(10,nhifr);

latex_hname={'$\\Delta E$','$\\dot{E}_{Pcp}$','$\\dot{E}_{Rn}$', ...
             '$\\dot{E}_{\\rho}$','$\\dot{E}_P$','$\\dot{H}+\\dot{L}$',...
             '$\\dot{E}_{DR}$','$\\dot{E}_{RO}$','$\\Delta C$','$\\dot{C}_{NEP}$'};
latex_hunit={'$GJ/m^2$','$GJ/m^2$','$GJ/m^2$','$GJ/m^2$', ...
              '$GJ/m^2$','$GJ/m^2$','$GJ/m^2$','$GJ/m^2$','$umol/m^2$','$umol/m^2$'};

for ih = 1:nhifr
    if(~(hpass(ih,2) && hpass(ih,3)))
        display(sprintf('Site: %s did not complete both main and test',...
                        hifrid{ih}));
    else
        display(sprintf('%s\n',hifrid{ih}));
        
        test_b_pfx = sprintf('%s/F_test_%s/test_%s_budget_state_patch_', ...
                                test_name,hifrid{ih},hifrid{ih});
        cont_b_pfx = sprintf('%s/F_main_%s/main_%s_budget_state_patch_', ...
                                test_name,hifrid{ih},hifrid{ih});

        display(test_b_pfx);

        id=strfind(test_b_pfx,'/');
        test_b_dir = test_b_pfx(1:id(end));
        
        id=strfind(cont_b_pfx,'/');
        cont_b_dir = cont_b_pfx(1:id(end));
        
        test_b_flist = dir(strcat(test_b_pfx,'*txt'));
        cont_b_flist = dir(strcat(cont_b_pfx,'*txt'));
        
        nbfiles     = length(test_b_flist);
        
        if (nbfiles ~= length(cont_b_flist))
            display(sprintf('Budget file lists are different lengths - %s',...
                            siteid{is}));
            display(['Different number of patch output, or multiple ' ...
                     'runs?']);
            return;
        elseif (nbfiles == 0)
            
            display('Could not find any patch files');
            return;
            
        else
            npatch = nbfiles;
            ndat1=0;
            for ipa=1:npatch
                tbfile = strcat(test_b_dir,test_b_flist(ipa).name);
                cbfile = strcat(cont_b_dir,cont_b_flist(ipa).name);
            
                % Read in the patch data for the test sim
                
                %[tdata_t,cdata_t,edata_t,wdata_t] =
                %read_patch_budgets(tbfile);
                
                [ndat,cbud_t,ebud_t,wbud_t,cstor_t,estor_t,wstor_t,~] ...
                    = read_patch_budgets(tbfile,ndat1);
                
                % Read in the patch data for the main sim
                [~,cbud_c,ebud_c,wbud_c,cstor_c,estor_c,wstor_c,dnv] ...
                    = read_patch_budgets(cbfile,ndat);
                
                % First pass, zero patch arrays
                if(ipa==1)
                    cbuds_t = zeros(ndat,5,npatch);
                    ebuds_t = zeros(ndat,9,npatch);
                    wbuds_t = zeros(ndat,7,npatch);
                    cstors_t = zeros(ndat,npatch);
                    estors_t = zeros(ndat,npatch);
                    wstors_t = zeros(ndat,npatch);
                    cbuds_c = zeros(ndat,5,npatch);
                    ebuds_c = zeros(ndat,9,npatch);
                    wbuds_c = zeros(ndat,7,npatch);
                    cstors_c = zeros(ndat,npatch);
                    estors_c = zeros(ndat,npatch);
                    wstors_c = zeros(ndat,npatch);
                end
                
                cbuds_t(:,:,ipa) = cbud_t;
                ebuds_t(:,:,ipa) = ebud_t;
                wbuds_t(:,:,ipa) = wbud_t;
                cbuds_c(:,:,ipa) = cbud_c;
                ebuds_c(:,:,ipa) = ebud_c;
                wbuds_c(:,:,ipa) = wbud_c;
                cstors_t(:,ipa) = cstor_t;
                estors_t(:,ipa) = estor_t;
                wstors_t(:,ipa) = wstor_t;
                cstors_c(:,ipa) = cstor_c;
                estors_c(:,ipa) = estor_c;
                wstors_c(:,ipa) = wstor_c;
                
                dtfac = (dnv(2)-dnv(1))*86400.0;
                ebud_cst=1e-6*dtfac*cumsum(ebud_t,1);
                cbud_cst=dtfac*cumsum(cbud_t,1);
                ebud_csc=1e-6*dtfac*cumsum(ebud_c,1);
                cbud_csc=dtfac*cumsum(cbud_c,1);
                
                latex_htab(1,ih) = latex_htab(1,ih)+...
                    (ebud_cst(end,2)-ebud_csc(end,2))./npatch;
                latex_htab(2,ih) = latex_htab(2,ih)+...
                    (ebud_cst(end,3)-ebud_csc(end,3))./npatch;
                latex_htab(3,ih) = latex_htab(3,ih)+...
                    (ebud_cst(end,4)-ebud_csc(end,4))./npatch;
                latex_htab(4,ih) = latex_htab(4,ih)+...
                    (ebud_cst(end,5)-ebud_csc(end,5))./npatch;
                latex_htab(5,ih) = latex_htab(5,ih)+...
                    (ebud_cst(end,6)-ebud_csc(end,6))./npatch;
                latex_htab(6,ih) = latex_htab(6,ih)+...
                    (ebud_cst(end,7)-ebud_csc(end,7))./npatch;
                latex_htab(7,ih) = latex_htab(7,ih)+...
                    (ebud_cst(end,8)-ebud_csc(end,8))./npatch;
                latex_htab(8,ih) = latex_htab(8,ih)+...
                    (ebud_cst(end,9)-ebud_csc(end,9))./npatch;
                latex_htab(9,ih) = latex_htab(9,ih)+...
                    (cbud_cst(end,2)-cbud_csc(end,2))./npatch;
                latex_htab(10,ih)= latex_htab(10,ih)+...
                    (cbud_cst(end,3)-cbud_csc(end,3))./npatch;
                
               
                ndat1=ndat;
            end
            
            
            
            % HIFI
            % INT_RESID INT_DEL_ES INT_PCP INT_RNET INT_RHO INT_P INT_ATMFL INT_DR
            % INT_RUN
            
            % ebuds_t = zeros(ndat,9,npatch);
            % EBUDS
            
            % Generate Plots
            % Time Series plot of patch 1 balance components
            
            cbudg_outfile{ih} = ...
                sprintf('%s/cbudg_%s_%s.eps',outdir,test_name,hifrid{ih});          
            ebudg_outfile{ih} = ...
                sprintf('%s/ebudg_%s_%s.eps',outdir,test_name,hifrid{ih});
            wbudg_outfile{ih} = ...
                sprintf('%s/wbudg_%s_%s.eps',outdir,test_name,hifrid{ih});
            
            
            plot_patchbudgets(dnv,cbuds_t,ebuds_t,wbuds_t,...
                                cstors_t,estors_t,wstors_t,...
                                cbuds_c,ebuds_c,wbuds_c,...
                                cstors_c,estors_c,wstors_c,...
                                cbudg_outfile{ih}, ...
                                ebudg_outfile{ih}, ...
                                wbudg_outfile{ih}, ...
                                visible);
            
            
            
            % Time Series plot of balance error
            
            % Summary Statistics
            
            
        end
    end        
end % for ih
end %if nhifr>0


%==========================================================================
% Validation of the gridded simulation
%==========================================================================

latex_gname = {'$AGB$','$BA$'};
latex_gunit = {'$kgC/ha$','$m^2/ha$'};


if ngrid>0
latex_gtab = zeros(2,ngrid);
   display('Assessing Gridded site(s)'); 
    for ig=1:ngrid
        
        if(~(gpass(ig,2) && gpass(ig,3)))
            display(sprintf('Grid: %s did not complete for main and test',...
                gridid{ig}));
        else
            
            display(sprintf('%s\n',gridid{ig}));

            test_gs_pfx = sprintf('%s/F_test_%s/test_%s-Q-',test_name,gridid{ig},gridid{ig});
            cont_gs_pfx = sprintf('%s/F_main_%s/main_%s-Q-',test_name,gridid{ig},gridid{ig});

            test_gs_flist = dir(strcat(test_gs_pfx,'*h5'));
            cont_gs_flist = dir(strcat(cont_gs_pfx,'*h5'));
            id=strfind(cont_gs_pfx,'/');
            cont_gs_dir = cont_gs_pfx(1:id(end));
            
            id=strfind(test_gs_pfx,'/');
            test_gs_dir = test_gs_pfx(1:id(end));
            
            tsfile = strcat(test_gs_dir,test_gs_flist(end).name);
            csfile = strcat(cont_gs_dir,cont_gs_flist(end).name);
            
            lat_gt    = double(hdf5read(tsfile,'/LATITUDE'));
            lon_gt    = double(hdf5read(tsfile,'/LONGITUDE'));
            pysiid_gt = hdf5read(tsfile,'/PYSI_ID');
            pysin_gt  = hdf5read(tsfile,'/PYSI_N');
            sipaid_gt = hdf5read(tsfile,'/SIPA_ID');
            sipan_gt  = hdf5read(tsfile,'/SIPA_N');
            pacoid_gt = hdf5read(tsfile,'/PACO_ID');
            pacon_gt  = hdf5read(tsfile,'/PACO_N');
            pft_gt    = double(hdf5read(tsfile,'/PFT'));
            areapa_gt = double(hdf5read(tsfile,'/AREA'));
            areasi_gt = double(hdf5read(tsfile,'/AREA_SI'));
            laico_gt  = double(hdf5read(tsfile,'/LAI_CO'));
            agbraw_gt = double(hdf5read(tsfile,'/AGB_PY')); %kg/m2 -> kgC/m2
            npoly_gt = length(lat_gt);
            
            lat_gc    = double(hdf5read(csfile,'/LATITUDE'));
            lon_gc    = double(hdf5read(csfile,'/LONGITUDE'));
            pysiid_gc = hdf5read(csfile,'/PYSI_ID');
            pysin_gc  = hdf5read(csfile,'/PYSI_N');
            sipaid_gc = hdf5read(csfile,'/SIPA_ID');
            sipan_gc  = hdf5read(csfile,'/SIPA_N');
            pacoid_gc = hdf5read(csfile,'/PACO_ID');
            pacon_gc  = hdf5read(csfile,'/PACO_N');
            pft_gc    = double(hdf5read(csfile,'/PFT'));
            areapa_gc = double(hdf5read(csfile,'/AREA'));
            areasi_gc = double(hdf5read(csfile,'/AREA_SI'));
            laico_gc = double(hdf5read(csfile,'/LAI_CO'));
            agbraw_gc = double(hdf5read(csfile,'/AGB_PY'));  %kg/m2 -> kgC/m2
            npoly_gc = length(lat_gc);
            
            if(npoly_gt~=npoly_gc)
                display('YOU SCREWED UP_ THE DATA');
                return;
            end
            
            [npft,ndbh,npoly] = size(agbraw_gt);
            
            lai_gt = zeros(npoly,npft);
            lai_gc = zeros(npoly,npft);
            agb_gt = zeros(npoly,npft);
            agb_gc = zeros(npoly,npft);
            
            for ipy=1:npoly
                isi_b = pysiid_gt(ipy);
                isi_e = isi_b + pysin_gt(ipy)-1;
                
                for isi=isi_b:isi_e
                    ipa_b = sipaid_gt(isi);
                    ipa_e = ipa_b + sipan_gt(isi)-1;
                    
                    for ipa=ipa_b:ipa_e
                        ico_b = pacoid_gt(ipa);
                        ico_e = ico_b + pacon_gt(ipa)-1;
                        for ico=ico_b:ico_e
                            ipft = pft_gt(ico);
                            lai_gt(ipy,ipft) = lai_gt(ipy,ipft)+...
                                laico_gt(ico)*areapa_gt(ipa)*areasi_gt(isi);
                        end
                    end
                end
                
                isi_b = pysiid_gc(ipy);
                isi_e = isi_b + pysin_gc(ipy)-1;
                for isi=isi_b:isi_e
                    ipa_b = sipaid_gc(isi);
                    ipa_e = ipa_b + sipan_gc(isi)-1;
                    for ipa=ipa_b:ipa_e
                        ico_b = pacoid_gc(ipa);
                        ico_e = ico_b + pacon_gc(ipa)-1;
                        for ico=ico_b:ico_e
                            ipft = pft_gc(ico);
                            lai_gc(ipy,ipft) = lai_gc(ipy,ipft)+...
                                laico_gc(ico)*areapa_gc(ipa)*areasi_gc(isi);
                        end
                    end
                end
                
                for ipft=1:npft
                    agb_gt(ipy,ipft) = sum(agbraw_gt(ipft,:,ipy));
                    agb_gc(ipy,ipft) = sum(agbraw_gc(ipft,:,ipy));
                end
                
            end
            
            % Determine which PFT's are present
            usepft = zeros(npft,1);
            for ipft=1:npft
                if(~isempty(find(agb_gt(:,ipft)>0)) || ...
                        ~isempty(find(agb_gc(:,ipft)))) %#ok<*EFIND>
                    usepft(ipft) = 1;
                end
            end
            
            agbmap_img{ig} = sprintf('%sagbmap_%s.eps',outdir,gridid{ig});
            laimap_img{ig} = sprintf('%slaimap_%s.eps',outdir,gridid{ig});
            
            plot_agbmaps(usepft,agb_gt,agb_gc,lon_gc,lat_gc,npoly,agbmap_img{ig},visible,grid_name{ig});
            plot_laimaps(usepft,lai_gt,lai_gc,lon_gc,lat_gc,npoly,laimap_img{ig},visible,grid_name{ig});
            
            latex_gtab(1,ig) = mean(sum(agb_gt,2)-sum(agb_gc,2));
            latex_gtab(2,ig) = mean(sum(lai_gt,2)-sum(lai_gc,2));
            
        end
    end   % for ig
end   % if ngrid>0


%==========================================================================
% Part 4: Table Fast Values
%==========================================================================



% HIFI
% INT_RESID INT_DEL_ES INT_PCP INT_RNET INT_RHO INT_P INT_ATMFL INT_DR
% INT_RUN

% GRID
% SUM_DEL_AGBT  SUM_DEL_AGB1 SUM_DEL_AGB2 ... SUM_DEL_AGBT







%==========================================================================
% Generate the Report
%==========================================================================

latexgen;
