myplot = gcf;
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;

name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

try

saveas(myplot,fullfile(pathstr, name),'fig')
catch
    error('err')
end