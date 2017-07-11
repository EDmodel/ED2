function [ndat,cbudget,ebudget,wbudget,cstorage,estorage,wstorage,dnv] ...
         =read_patch_budgets(filename,in_nline)
    
%=====================================================================
% Read in and parse EDM patch budget data
%=====================================================================
    
    
% Text input

% 1YEAR 2MONTH 3DAY 4TIME 5LAI 6WAI 7HEIGHT 8CO2.STORAGE 9CO2.RESIDUAL 
% 10CO2.DSTORAGE 11CO2.NEP 12CO2.DENS.EFF 13CO2.LOSS2ATM 14ENE.STORAGE 
% 15ENE.RESIDUAL 16ENE.DSTORAGE 17ENE.PRECIP 18ENE.NETRAD 19ENE.DENS.EFF 
% 20ENE.PRSS.EFF 21ENE.LOSS2ATM 22ENE.DRAINAGE 23ENE.RUNOFF 24H2O.STORAGE 
% 25H2O.RESIDUAL 26H2O.DSTORAGE 27H2O.PRECIP 28H2O.DENS.EFF
% 29H2O.LOSS2ATM 30H2O.DRAINAGE 31H2O.RUNOFF
    
% First lets get the number of lines in the file
    fid=fopen(filename);
    ndat=0;
    fgetl(fid);
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        ndat=ndat+1;
    end
    fclose(fid);
       
    % Now lest initialize the output arrays
    % and read in the data
    
    if(in_nline ~= 0 && ndat~=in_nline)
        display(['Number of entries in patch budgets does not ' ...
                 'match']);
        return;
    else

        % Initialize output arrays
        
        %        yearv  = zeros(ndat,1); %1
        %        monthv = zeros(ndat,1); %2
        %        dayv   = zeros(ndat,1); %3
        %        timev  = zeros(ndat,1); %4
        dnv    = zeros(ndat,1);
        
        laiv   = zeros(ndat,1); %5
        waiv   = zeros(ndat,1); %6
        hgtv   = zeros(ndat,1); %7
        
        cstorage= zeros(ndat,1);
        cbudget = zeros(ndat,5);
        estorage= zeros(ndat,1);
        ebudget = zeros(ndat,9);
        wstorage= zeros(ndat,1);
        wbudget = zeros(ndat,7);
        
 %       cstor  = zeros(ndat,1); %8
 %       cres   = zeros(ndat,1); %9
 %       cdstor = zeros(ndat,1); %10
 %       cnep   = zeros(ndat,1); %11
 %       cdens  = zeros(ndat,1); %12
 %       catm   = zeros(ndat,1); %13

 %       estor  = zeros(ndat,1); %14
 %       eres   = zeros(ndat,1); %15
 %       edstor = zeros(ndat,1); %16
 %       epcp   = zeros(ndat,1); %17
 %       ernet  = zeros(ndat,1); %18
 %       edens  = zeros(ndat,1); %19
 %       eprss  = zeros(ndat,1); %20
 %       eatm   = zeros(ndat,1); %21
 %       edrain = zeros(ndat,1); %22
 %       erunn  = zeros(ndat,1); %23
        
 %       wstor  = zeros(ndat,1); %24
 %       wres   = zeros(ndat,1); %25
 %       wdstor = zeros(ndat,1); %26
 %       wpcp   = zeros(ndat,1); %27
 %       wdens  = zeros(ndat,1); %28
 %       watm   = zeros(ndat,1); %29
 %       wdrain = zeros(ndat,1); %30
 %       wrunn  = zeros(ndat,1); %31
        
        fid=fopen(filename);
        fgetl(fid);
        id=0;
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            id=id+1;
            dvec = str2num(tline);
            
            yearv  = dvec(1);
            monthv = dvec(2);
            dayv   = dvec(3);
            timev  = dvec(4);
            
            dnv(id) = datenum(yearv,monthv,dayv,0,0,timev);
        
            laiv(id)   = dvec(5);
            waiv(id)   = dvec(6); %6
            hgtv(id)   = dvec(7); %7
            cstorage(id)     = dvec(8);
            cbudget(id,1:5)  = dvec(9:13);
            estorage(id)     = dvec(14);
            ebudget(id,1:9) = dvec(15:23);
            wstorage(id)     = dvec(24);
            wbudget(id,1:7)  = dvec(25:31);
            
           
            
        end
        fclose(fid);
        
    end
    

end
