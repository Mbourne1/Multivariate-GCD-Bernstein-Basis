global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
       
       
        % Plot 3d data points
        figure_name = sprintf('%s - 3d plot',mfilename);
        figure('name',figure_name)
        hold on
        title('Minimum Singular Values in S_{t_{1},t_{2}}')
        xlabel('t_{1}')
        ylabel('t_{2}')
        scatter3(data(:,1),data(:,2),log10(data(:,4)),'filled')
        grid('on')
        hold off
        
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end