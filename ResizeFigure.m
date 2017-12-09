
myplot = gca;

myval_side = 0.10;
myval_base = 0.08;

% Use for plotting 3d graphics
%
% xlabel('k_{1}')
% ylabel('k_{2}')
% zlabel('log_{10} (\sigma_{k_{1}, k_{2}})')


% xlabel('j : column index')
% ylabel('Coefficient Multiplier')


%ylabel('$\log_{10}\left( \Re \right)$')
%xlabel('$i$')

% Use for singular value plotting
%
%xlabel('$k$')
%ylabel('$\log_{10}\left( \sigma_{k,i} \right)$')

% Use for minimum singular values
%xlabel('$k$','Interpreter','LaTex')
%ylabel('$\log_{10} \left( \sigma_{k} \right)$','Interpreter','LaTex')


% Label for DTQ
%ylabel('$\log_{10}\left( {{m}\choose{i}} {{n-k}\choose{j}} / {{m + n - k}\choose{i + j}} \right)$', 'Interpreter','LaTex')

% Label for DT
%ylabel('$\log_{10}\left( {{m}\choose{i}}  / {{m + n - k}\choose{i + j}} \right)$', 'Interpreter','LaTex')

% Label for DTQ Denom removed
%ylabel('$\log_{10}\left( {{i + j}\choose{i}} {{m + n - k -i -j}\choose{m - i}} \right)$', 'Interpreter','LaTex')

%xlabel('$j : Column Index$','Interpreter','LaTex')

% Label for TQ
%ylabel('$\log_{10}\left( {{m}\choose{i}} {{n-k}\choose{j}} \right)$', 'Interpreter','LaTex')

% Use for residuals in stln
%xlabel('i : Iteration')
%ylabel('log_{10} (Residual)')


hLines = findobj(gca,'Type','line');



for k = 1:length(hLines)
    activeLine = hLines(k);
    
    %set(activeLine,'LineWidth',2)
    %set(activeLine,'markers',4)
    %set(activeLine,'Marker','*')
    %set(activeLine,'Marker','square')
    %set(activeLine, 'MarkerEdgeColor','black')
    %set(activeLine,'MarkerFaceColor','black')
    
end


%set(gca,'FontSize',20)
set(gca,'FontName','Helvetica')

set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])


myplot = gcf;
grid on
box on
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;




name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

%print(gcf,'-depsc',fullfile(pathstr, name1));



name2 = strcat(name,'.fig');
saveas(myplot,fullfile(pathstr, name2),'fig')


%name2 = strcat(name,'.png');
%saveas(myplot,fullfile(pathstr, name2),'png')