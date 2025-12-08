%
% script called by dev_test to generate reports
% repgen.m
%
%

% Create RptgenML.CReport
RptgenML_CReport1 = RptgenML.CReport('Description',...
    'EDM 2.2 automated verification of a new model commit.',...
    'Format','pdf-fop','Stylesheet','simple_landscape');


% Create the title page

report_title_str = sprintf(...
    'EDM 2.2 automated verification of commit: %s',test_name);
if(strcmp(tester_name,committer_name))
    authors = tester_name;
else
    authors = sprintf('%s and %s',tester_name,committer_name);
end
rptgen_cfr_titlepage1 = rptgen.cfr_titlepage('AuthorMode','manual',...
    'Author',authors,'Title',report_title_str,...
    'DoSinglePage',true);
setParent(rptgen_cfr_titlepage1,RptgenML_CReport1);


% Create the test specifications page

spec_text = sprintf('Test Version Branched from: %s\n', ...
                           branch_version);
spec_text = sprintf('%sCommitter (changed model code): %s\n', ...
                     spec_text,committer_name);
spec_text = sprintf('%sTester (generated this report): %s\n\n', ...
                     spec_text,tester_name);
spec_text = sprintf('%sDescription of Changes: %s',...
                     spec_text,test_description);

rptgen_cfr_section0 = rptgen.cfr_section('SectionTitle',...
    'Test Specification');
rptgen_cfr_textpage0 = rptgen.cfr_text('isLiteral',true,...
    'Content', spec_text);
setParent(rptgen_cfr_textpage0,rptgen_cfr_section0);
setParent(rptgen_cfr_section0,RptgenML_CReport1);

% Create the verifications page

rptgen_cfr_sectionA = rptgen.cfr_section('SectionTitle',...
    'Verification of Completion');
setParent(rptgen_cfr_sectionA,RptgenML_CReport1);

verifications = sprintf('\n');
verifications = [verifications,sprintf('The following simulations were executed in debug-mode\n\n')];

for is=1:nsite
    if(spass(is,1)==1)
        verifications = [verifications,sprintf('%s: COMPLETED\n',siteid{is})]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('%s: DID NO COMPLETE\n',siteid{is})]; %#ok<AGROW>
    end
end

for ig=1:ngrid
    if(gpass(ig,1)==1)
        verifications = [verifications,sprintf('%s: COMPLETED\n',gridid{ig})]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('%s: DID NOT COMPLETE\n',gridid{ig})]; %#ok<AGROW>
    end
end


rptgen_cfr_paragraphA = rptgen.cfr_paragraph;
rptgen_cfr_textA = rptgen.cfr_text('isLiteral',true,...
    'Content', verifications);
set(rptgen_cfr_paragraphA,'ParaTextComp',rptgen_cfr_textA);
setParent(rptgen_cfr_paragraphA,rptgen_cfr_sectionA);


% Create the Simulation Specifications

rptgen_cfr_sectionB = rptgen.cfr_section('SectionTitle',...
    'Simulation Specifications');
setParent(rptgen_cfr_sectionB,RptgenML_CReport1);

% Create that table that shows all the specs

rptgen_cfr_desc_img = rptgen.cfr_image('FileName',...
            'test_suit_combos_v4_cropped.png','DocHorizAlign','center');

setParent(rptgen_cfr_desc_img,rptgen_cfr_sectionB);


if nsite>0
    for is=1:nsite
        
        % First SOI
        % =========================================================================
        % Create rptgen.cfr_section
        rptgen_cfr_section(is) = rptgen.cfr_section('SectionTitle',site_name{is});
        setParent(rptgen_cfr_section(is),RptgenML_CReport1);
        
        if(~(spass(is,2) && spass(is,3)))
            
            rptgen_cfr_paragraph = rptgen.cfr_paragraph;
            rptgen_cfr_text = rptgen.cfr_text('isLiteral',true,'Content', ...
                sprintf('\n\n**** AT LEAST ONE SIMULATION DID NOT COMPLETE \n\n'));
            set(rptgen_cfr_paragraph,'ParaTextComp',rptgen_cfr_text);
            setParent(rptgen_cfr_paragraph,rptgen_cfr_section(is));
            
        else
            
            % Subsection Fluxes
            rptgen_cfr_section1(is) = rptgen.cfr_section('SectionTitle','Biophysical 1');
            setParent(rptgen_cfr_section1(is),rptgen_cfr_section(is));
            
            % Create rptgen.cfr_image
            rptgen_cfr_fluximage(is) = rptgen.cfr_image('FileName',fluxes_img{is},...
                                       'DocHorizAlign','center');
            setParent(rptgen_cfr_fluximage(is),rptgen_cfr_section1(is));
            
            % Subsection Fluxes
            rptgen_cfr_section12(is) = rptgen.cfr_section('SectionTitle','Biophysical 2');
            setParent(rptgen_cfr_section12(is),rptgen_cfr_section(is));
            
            % Create rptgen.cfr_image
            rptgen_cfr_stateimage(is) = rptgen.cfr_image('FileName',states_img{is},...
                                       'DocHorizAlign','center');
            setParent(rptgen_cfr_stateimage(is),rptgen_cfr_section12(is));
             
            % Subsection Succession
            if(plongterm_plt(is)>0)
                rptgen_cfr_section2(is) = rptgen.cfr_section('SectionTitle','Succession');
                setParent(rptgen_cfr_section2(is),rptgen_cfr_section(is));
                
                % Create rptgen.cfr_image
                rptgen_cfr_ltermimage(is) = rptgen.cfr_image('FileName',longterm_img{is},...
                    'DocHorizAlign','center'); %...
 %                   'ViewPortType','fixed',...
 %                   'ViewPortSize',[8*1.4 6*1.4]);
                setParent(rptgen_cfr_ltermimage(is),rptgen_cfr_section2(is));
                
            end
            
            % Subsection Ecosystem Profiles
            rptgen_cfr_section3(is) = rptgen.cfr_section('SectionTitle','End Simulation Ecosystem Profiles');
            setParent(rptgen_cfr_section3(is),rptgen_cfr_section(is));
            
            % Create rptgen.cfr_image
            rptgen_cfr_profimage(is) = rptgen.cfr_image('FileName',strcomp_img{is},...
                'DocHorizAlign','center'); %#ok<*SAGROW>
            setParent(rptgen_cfr_profimage(is),rptgen_cfr_section3(is));
            
        end
        
    end % for nsite
end % if nsite>0

% Hi Frequency Output (Energy,Carbon and Mass Balances)
% =========================================================================

if nhifr>0
    for ih=1:nhifr
        
        isec = ih+nsite;
        
        % First HIFR
        % =========================================================================
        % Create rptgen.cfr_section
        rptgen_cfr_section(isec) = rptgen.cfr_section('SectionTitle',hifr_name{ih});
        setParent(rptgen_cfr_section(isec),RptgenML_CReport1);
        
        if(~(hpass(ih,2) && hpass(ih,3)))
            
            rptgen_cfr_paragraph = rptgen.cfr_paragraph;
            rptgen_cfr_text = rptgen.cfr_text('isLiteral',true,'Content', ...
                sprintf('\n\n**** AT LEAST ONE SIMULATION DID NOT COMPLETE \n\n'));
            set(rptgen_cfr_paragraph,'ParaTextComp',rptgen_cfr_text);
            setParent(rptgen_cfr_paragraph,rptgen_cfr_section(isec));
            
        else
            
            % Subsection Energy Budget
            rptgen_cfr_section1(isec) = rptgen.cfr_section('SectionTitle',...
                'Energy Budget');
            setParent(rptgen_cfr_section1(isec),rptgen_cfr_section(isec));
            
            % Create Image
            rptgen_cfr_ebudgimage(isec) = rptgen.cfr_image('FileName',...
                ebudg_outfile{ih},'DocHorizAlign','center');
            setParent(rptgen_cfr_ebudgimage(isec),rptgen_cfr_section1(isec));
            
            % Subsection Water Mass Budget
            rptgen_cfr_section2(isec) = rptgen.cfr_section('SectionTitle',...
                'Water Mass Budget');
            setParent(rptgen_cfr_section2(isec),rptgen_cfr_section(isec));
            
            % Create Image
            rptgen_cfr_wbudgimage(isec) = rptgen.cfr_image('FileName',...
                wbudg_outfile{ih},'DocHorizAlign','center');
            setParent(rptgen_cfr_wbudgimage(isec),rptgen_cfr_section2(isec));

            % Subsection Carbon Budget
            rptgen_cfr_section3(isec) = rptgen.cfr_section('SectionTitle',...
                'Carbon Budget');
            setParent(rptgen_cfr_section3(isec),rptgen_cfr_section(isec));
            
            % Create Image
            rptgen_cfr_cbudgimage(isec) = rptgen.cfr_image('FileName',...
                cbudg_outfile{ih},'DocHorizAlign','center');
            setParent(rptgen_cfr_cbudgimage(isec),rptgen_cfr_section3(isec));
        end
    end
        
end
    
    
% Offline Gridded Run
% =========================================================================
% Create rptgen.cfr_section

if ngrid>0
    for ig=1:ngrid
        
        rptgen_cfr_section_g = rptgen.cfr_section('SectionTitle',grid_name{ig});
        setParent(rptgen_cfr_section_g,RptgenML_CReport1);
        
        if(~(gpass(ig,2) && gpass(ig,3)))
            
            rptgen_cfr_paragraph = rptgen.cfr_paragraph;
            rptgen_cfr_text = rptgen.cfr_text('isLiteral',true,'Content', ...
                sprintf('\n\n**** AT LEAST ONE SIMULATION DID NOT COMPLETE \n\n'));
            set(rptgen_cfr_paragraph,'ParaTextComp',rptgen_cfr_text);
            setParent(rptgen_cfr_paragraph,rptgen_cfr_section_g);
            
        else
            
            % Subsection Fluxes
            rptgen_cfr_section_g1 = rptgen.cfr_section('SectionTitle','AGB Maps');
            setParent(rptgen_cfr_section_g1,rptgen_cfr_section_g);
            
            
            % Create rptgen.cfr_image
            rptgen_cfr_agbimage = rptgen.cfr_image('FileName',agbmap_img{ig},...
                'DocHorizAlign','center');
            
            setParent(rptgen_cfr_agbimage,rptgen_cfr_section_g1);
            
            % Subsection Fluxes
            rptgen_cfr_section_g2 = rptgen.cfr_section('SectionTitle','LAI Maps');
            setParent(rptgen_cfr_section_g2,rptgen_cfr_section_g);
            
            % Create rptgen.cfr_image
            rptgen_cfr_laiimage = rptgen.cfr_image('FileName',laimap_img{ig},...
                'DocHorizAlign','center'); 
            setParent(rptgen_cfr_laiimage,rptgen_cfr_section_g2);
            
        end
    end %for ig
end %if ngrid>0

pdf_file_name = strcat('-o',outdir,test_name,'.pdf');
[report1]=report(RptgenML_CReport1,pdf_file_name,'-noview','-fPDF')

