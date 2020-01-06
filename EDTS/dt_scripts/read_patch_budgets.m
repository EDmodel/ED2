function [ndat,cbudget,ebudget,wbudget,cstorage,estorage,wstorage,dnv] ...
         =read_patch_budgets(filename,in_nline)
    
%=====================================================================
% Read in and parse EDM patch budget data
%=====================================================================
    
    
% Text input
%  current_time%year      , current_time%month     , current_time%date       3
%  current_time%time      , patch_lai              , patch_wai               6
%  csite%veg_height(ipa)  , co2budget_finalstorage , co2curr_residual        9
%  co2budget_deltastorage , co2curr_nep            , co2curr_denseffect     12
%  co2curr_zcaneffect     , co2curr_loss2atm       , cbudget_finalstorage   15
%  cbudget_committed      , ccurr_residual         , cbudget_deltastorage   18
%  ccurr_denseffect       , ccurr_zcaneffect       , ccurr_seedrain         21
%  ccurr_loss2yield       , ccurr_loss2atm         , ebudget_finalstorage   24
%  ecurr_residual         , ebudget_deltastorage   , ecurr_precipgain       27
%  ecurr_netrad           , ecurr_denseffect       , ecurr_prsseffect       30
%  ecurr_hcapeffect       , ecurr_wcapeffect       , ecurr_zcaneffect       33
%  ecurr_pheneffect       , ecurr_loss2atm         , ecurr_loss2drainage    36
%  ecurr_loss2runoff      , wbudget_finalstorage   , wcurr_residual         39
%  wbudget_deltastorage   , wcurr_precipgain       , wcurr_denseffect       42
%  wcurr_wcapeffect       , wcurr_zcaneffect       , wcurr_pheneffect       45
%  wcurr_loss2atm         , wcurr_loss2drainage    , wcurr_loss2runoff      48
    
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
        cbudget = zeros(ndat,6);
        estorage= zeros(ndat,1);
        ebudget = zeros(ndat,9);
        wstorage= zeros(ndat,1);
        wbudget = zeros(ndat,7);
        
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
        
            laiv(id)          = dvec(5);
            waiv(id)          = dvec(6); %6
            hgtv(id)          = dvec(7); %7
            % Standardise carbon output.
            cstorage(id)      = dvec(15)+dvec(16);
            cbudget(id,1:4)   = dvec(17:20);
            cbudget(id,5)     = dvec(21)-dvec(22);
            cbudget(id,6)     = -dvec(23)
            % Simplify energy output.
            estorage(id)      = dvec(24);
            ebudget(id,1:6)   = dvec(25:30);
            ebudget(id,7)     = dvec(31)+dvec(32)+dvec(33)+dvec(34);
            ebudget(id,8)     = -dvec(35);
            ebudget(id,9)     = dvec(36)+dvec(37);
            % Simplify water output.
            wstorage(id)      = dvec(38);
            wbudget(id,1:4)   = dvec(39:42);
            wbudget(id,5)     = dvec(43)+dvec(44)+dvec(45);
            wbudget(id,6)     = -dvec(46);
            wbudget(id,7)     = dvec(47)+dvec(48);

        end
        fclose(fid);
        
    end
    

end
