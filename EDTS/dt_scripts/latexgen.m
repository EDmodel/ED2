%==========================================================================
%
% Generate a latex beamer report from matlab output
%
%
%
%==========================================================================



report_filename = strcat(outdir,test_name);

report_texname  = strcat(report_filename,'.tex');
report_dviname  = strcat(report_filename,'.dvi');
report_pdfname  = strcat(report_filename,'.pdf');

% Start off by generating the header information from the template

display(report_texname);

cmdstr = sprintf('cp %s %s', ...
    strcat(pwd,'/dt_scripts/latextemplates/template_p1.tex'), ...
    report_texname);

system(cmdstr);

% Copy the beaver
cmdstr = sprintf('cp %s %s', ...
    strcat(pwd,'/dt_scripts/latextemplates/beamercolorthemebeaver.sty'), ...
    strcat(outdir,'/'));

system(cmdstr);

% Write out bunch of the header stuff

fid = fopen(report_texname,'a');

% Title
test_nametex=strrep(test_name,'_','$\_$')
textinsrt = sprintf('EDM 2.2 automated verification of commit: %s',test_nametex);
textinsrt = sprintf('title[%s]{%s}',test_nametex,textinsrt);
fprintf(fid,'\\%s\n',textinsrt);


%Author
if(strcmp(tester_name,committer_name))
    authors = tester_name;
else
    authors = sprintf('%s and %s',tester_name,committer_name);
end
textinsrt = sprintf('author{%s}',authors);
fprintf(fid,'\\%s\n',textinsrt);


%Date
%textinsrt = sprintf('date');
fprintf(fid,'\\date{\\today}\n');


% Begin document
fprintf(fid,'\\begin{document}\n');

% Frame the title-page
fprintf(fid,'\\frame{\\titlepage}\n');

% Frame the table of contents
fprintf(fid,'\\frame{\\frametitle{Table of Contents}\\scriptsize\\tableofcontents}\n');




% Section 1 - Testing Specifications
% =========================================================================

fprintf(fid,'\\section{Test Specifications}\n');

fprintf(fid,'\\scriptsize\n');

textinsrt = sprintf('Test Version Branched from: %s\\\\ \n',branch_version);
textinsrt = sprintf('%sCommitter (changed model code): %s\\\\ \n', ...
    textinsrt,committer_name);
textinsrt = sprintf('%sTester (generated this report): %s\\\\[0.5cm] \n', ...
    textinsrt,tester_name);

% Reformat latex special characters
textdesc=latex_prep(test_description);
textinsrt = sprintf('%sDescription of Changes: %s',...
    textinsrt,textdesc);
fprintf(fid,'\\frame{\n%s\n}\n',textinsrt);

fclose(fid);

% Generate the parameter evaluation tables

cmdstr = sprintf('./dt_scripts/gen_param_table.sh %s',test_name);
system(cmdstr);



% Concatenate the spec table 1
table1_texname = strcat(outdir,'/vartable1.tex');
cmdstr = sprintf('cat %s >> %s',table1_texname,report_texname);
system(cmdstr);

% Concatenate the spec table 2
table1_texname = strcat(outdir,'/vartable2.tex');
cmdstr = sprintf('cat %s >> %s',table1_texname,report_texname);
system(cmdstr);

% Concatenate the spec table 3
table1_texname = strcat(outdir,'/vartable3.tex');
cmdstr = sprintf('cat %s >> %s',table1_texname,report_texname);
system(cmdstr);

% Concatenate the spec table 4
table1_texname = strcat(outdir,'/vartable4.tex');
cmdstr = sprintf('cat %s >> %s',table1_texname,report_texname);
system(cmdstr);


fid = fopen(report_texname,'a+');

%\\includegraphics[width=0.99\\textwidth]{dt_scripts/test_suit_combos_v4_cropped.eps}}\n');


% Section 2 - Testing Summary
% =========================================================================

fprintf(fid,'\\section{Test Summary}\n');
fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\subsection{Debug Completion Test}\n');
verifications = sprintf('\n');
verifications = [verifications,sprintf('\\begin{scriptsize}\n')];
verifications = [verifications,sprintf('The following simulations resulted in completion or failure: \\\\[0.5cm] \n\n')];

verifications = [verifications,sprintf('RUN~ \\quad DBUG \\quad TEST \\quad  MAIN \\\\ \n\n')];

for is=1:nsite
    if(spass(is,1)==1)
        verifications = [verifications,sprintf('%s: \\quad COMP ',siteid{is})]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('%s: \\quad FAIL ',siteid{is})]; %#ok<AGROW>
    end
    if(spass(is,2)==1)
        verifications = [verifications,sprintf('\\quad COMP ')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL ')]; %#ok<AGROW>
    end
    if(spass(is,3)==1)
        verifications = [verifications,sprintf('\\quad COMP \\\\ \n')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL \\\\ \n')]; %#ok<AGROW>
    end

end

for ih=1:nhifr
    if(hpass(ih,1)==1)
        verifications = [verifications,sprintf('%s: \\quad COMP ',hifrid{ih})]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('%s: \\quad FAIL ',hifrid{ih})]; %#ok<AGROW>
    end
    if(hpass(ih,2)==1)
        verifications = [verifications,sprintf('\\quad COMP ')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL ')]; %#ok<AGROW>
    end
    if(hpass(ih,3)==1)
        verifications = [verifications,sprintf('\\quad COMP \\\\ \n')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL \\\\ \n')]; %#ok<AGROW>
    end


end

for ig=1:ngrid
    if(gpass(ig,1)==1)
        verifications = [verifications,sprintf('%s: \\quad COMP ',gridid{ig})]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('%s: \\quad FAIL ',gridid{ig})]; %#ok<AGROW>
    end
    if(gpass(ig,2)==1)
        verifications = [verifications,sprintf('\\quad COMP ')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL ')]; %#ok<AGROW>
    end
    if(gpass(ig,3)==1)
        verifications = [verifications,sprintf('\\quad COMP \\\\ \n')]; %#ok<AGROW>
    else
        verifications = [verifications,sprintf('\\quad FAIL \\\\ \n')]; %#ok<AGROW>
    end


end
verifications = [verifications,sprintf('\\end{scriptsize}\n')];


fprintf(fid,'\\frame{\n%s\n}\n',verifications);

fprintf(fid,'\\subsection{Output Comparison Table}\n');


% Test Summary - Table
fprintf(fid,'\\frame{\n');
fprintf(fid,'\\tiny\n');
fprintf(fid,'\\begin{columns}\n');
fprintf(fid,'\\begin{column}{\\textwidth}\n');
fprintf(fid,'\\vbox to .8\\textheight{ \n');
%fprintf(fid,'\\begin{table}\n');
fprintf(fid,'\\begin{center}\n');
    
if nsite>0
    fprintf(fid,'\\begin{tabular}{c} SOI Run(s) \\end{tabular}\n');
    fprintf(fid,'\\begin{tabular}{|r|c|c|c|c|c|c|c|c|c|}\n');
    fprintf(fid,'\\hline\n');
    textinsrt='Site ';
    for k=1:9
        textinsrt=[textinsrt,sprintf('& %s',latex_fname{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    textinsrt=' ';
    for k=1:9
        textinsrt=[textinsrt,sprintf('& %s',latex_funit{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    fprintf(fid,'\\hline\n');
    for is=1:nsite
	     if(spass(is,2) && spass(is,3))
        fprintf(fid,'%s & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n',...
            siteid{is},latex_ftab(1,is),latex_ftab(2,is),latex_ftab(3,is), ...
                       latex_ftab(3,is),latex_ftab(4,is),latex_ftab(5,is), ...
                       latex_ftab(6,is),latex_ftab(8,is),latex_ftab(9,is));
end
    end
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\vfill \n');
end



if nhifr>0
    fprintf(fid,'\\begin{tabular}{c} Hi-Frequency Run(s) (Time-Integrated) \\end{tabular}\n');
    fprintf(fid,'\\begin{tabular}{ |r|c|c|c|c|c|c|c|c|c|c|}\n');
    fprintf(fid,'\\hline \n');
    textinsrt='Site ';
    for k=1:10
        textinsrt=[textinsrt,sprintf('& %s',latex_hname{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    textinsrt=' ';
    for k=1:10
        textinsrt=[textinsrt,sprintf('& %s',latex_hunit{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    fprintf(fid,'\\hline\n');
    for ih=1:nhifr
		if(hpass(ih,2) && hpass(ih,3))
        fprintf(fid,'%s & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f\\\\ \n',...
            hifrid{ih},latex_htab(1,ih),latex_htab(2,ih),latex_htab(3,ih), ...
                       latex_htab(3,ih),latex_htab(4,ih),latex_htab(5,ih), ...
                       latex_htab(6,ih),latex_htab(8,ih),latex_htab(9,ih), ...
                       latex_htab(10,ih));
    end
end
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\vfill \n');
end

if ngrid>0
    fprintf(fid,'\\begin{tabular}{c} Gridded Run(s) \\end{tabular}\n');
    fprintf(fid,'\\\\');
    fprintf(fid,'\\begin{tabular}{ |r|c|c|}\n');
    fprintf(fid,'\\hline \n');
    textinsrt='Site ';
    for k=1:2
        textinsrt=[textinsrt,sprintf('& %s',latex_gname{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    textinsrt=' ';
    for k=1:2
        textinsrt=[textinsrt,sprintf('& %s',latex_gunit{k})];
    end
    textinsrt=[textinsrt,'\\\\ \n'];
    fprintf(fid,textinsrt);
    fprintf(fid,'\\hline\n');
    for ig=1:ngrid
	     if(gpass(ig,2) && gpass(ig,3))
        fprintf(fid,'%s & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f\\\\ \n',...
            gridid{ig},latex_gtab(1,ig),latex_gtab(2,ig));
end
    end
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\vfill \n');
end


%fprintf(fid,'\\caption{Simulation mean differences (TEST-MAIN).}\n');
fprintf(fid,'\\end{center}\n');
%fprintf(fid,'\\end{table}\n');
fprintf(fid,'} \n');
fprintf(fid,'\\end{column}\n');
fprintf(fid,'\\end{columns}\n');
fprintf(fid,'}\n');


% Sections 3+ SOIs
% =========================================================================

if nsite>0
    fprintf(fid,'\\section{\nSite of Interest (SOI) Runs\n}\n');
    for is=1:nsite
        fprintf(fid,'\\subsection{%s}\n',site_name{is});
        if(~(spass(is,2) && spass(is,3)))
            fprintf(fid,'\\frame{\n AT LEAST ONE SIMULATION DID NOT COMPLETE \\\\ COMPARATIVE ANALYSIS IMPOSSIBLE\n}\n');
        else
            
            % Ecosystem Profiles
            fprintf(fid,'\\frame{\\noindent\\includegraphics[width=0.48\\textwidth]{%s}\\includegraphics[width=0.48\\textwidth]{%s}}\n',strcomp_timg{is},strcomp_cimg{is});

            % Biophysics Variables 1
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',fluxes_img{is});
            
            % Biophysics Variables 2
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',states_img{is});
            
            % Succession
            fprintf(fid,'\\frame{\\noindent\\includegraphics[width=0.64\\textwidth]{%s}\\includegraphics[width=0.32\\textwidth]{%s}}\n',pftsucc_img{is},soilcarb_img{is});
        end
    end
end


% Sections for High Freuqency RUns
% =========================================================================

if nhifr>0
    fprintf(fid,'\\section{\nHigh Frequency Output\n}\n');
    for ih=1:nhifr
        fprintf(fid,'\\subsection{%s}\n',hifr_name{ih});
        if(~(hpass(ih,2) && hpass(ih,3)))
            fprintf(fid,'\\frame{\n AT LEAST ONE SIMULATION DID NOT COMPLETE \\\\ COMPARATIVE ANALYSIS IMPOSSIBLE\n}\n');
        else

            % Energy Budget
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',ebudg_outfile{ih});
            
            % Water Budget
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',wbudg_outfile{ih});
            
            % Carbon Budget
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',cbudg_outfile{ih});
            
        end
    end
end

% Sections for Gridded Runs
% =========================================================================

if ngrid>0
    fprintf(fid,'\\section{\nGridded Output\n}\n');
    for ig=1:ngrid
        fprintf(fid,'\\subsection{%s}\n',grid_name{ig});
        if(~(gpass(ig,2) && gpass(ig,3)))
            fprintf(fid,'\\frame{\n AT LEAST ONE SIMULATION DID NOT COMPLETE \\\\ COMPARATIVE ANALYSIS IMPOSSIBLE\n}\n');
        else
            % AGB Maps
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',agbmap_img{ig});
            
            % LAI Maps
            fprintf(fid,'\\frame{\\includegraphics[height=0.92\\textheight]{%s}}\n',laimap_img{ig});
            
        end
    end
end






% End the Document
textinsrt = sprintf('end{document}');
fprintf(fid,'\\%s',textinsrt);

fclose(fid);




% Compile the latex document
cmdstr = sprintf('latex --output-directory=%s  %s',outdir,report_texname);
system(cmdstr);
system(cmdstr); %Again for outline purposes

% Generate the pdf from the dvi

cmdstr = sprintf('dvipdf %s %s',report_dviname,report_pdfname);
system(cmdstr);



