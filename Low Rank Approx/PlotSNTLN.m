switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','STLN - Residuals in STLN')
        hold on
        title('Residuals in STLN')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
        plot((1:1:ite),log10(condition),'-s')
        hold off
        
       
    case 'n'
    otherwise
        error('err')
end