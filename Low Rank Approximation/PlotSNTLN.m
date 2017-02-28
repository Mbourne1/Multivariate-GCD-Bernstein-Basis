if(SETTINGS.PLOT_GRAPHS)
    
        figure('name','STLN - Residuals in STLN')
        hold on
        title('Residuals in STLN')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
        plot((1:1:ite),log10(condition),'-s')
        hold off
        
       
    
end