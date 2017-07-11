
% Some characters need the escape character in latex, this adds it before
% these various trouble characters

function outstr=latex_prep(instr)


trouble = {'#','$','%','&','\','^','_','{','}','~'};

fixer = {'\#','\$','\%','\&','\textbackslash{}',...
    '\textasciicircum{}','\_','\}','\textasciitilde{}'};

nt = length(trouble);

for it=1:nt

ids=regexp(instr,trouble{it});
n  = length(ids);

if(n>0)
    outstr = '';
    id1=1;
    for k=1:n
        id2=ids(k)-1;
        outstr = strcat(outstr,instr(id1:id2),fixer{it});
        id1=id2+2;
    end
    outstr=strcat(outstr,instr(id1:end));
else
    outstr=instr;
end

instr=outstr;

end