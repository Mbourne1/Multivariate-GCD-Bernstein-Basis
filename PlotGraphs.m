global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        %% plot the minimum singular values
        figure_name = sprintf('%s : Minimum Singular Values',mfilename);
        figure('name',figure_name)
        title('minimum Singular Value for each subresultant matrix S_{k,k}')
        hold on
        plot(log10(vMinimumSingularValues),'-s','DisplayName','Preprocessed');
        %plot(log10(min_sing_val_vec_unproc),'-s','DisplayName','Unprocessed');
        xlabel('k : index of subresultant')
        legend(gca,'show')
        ylabel('log_{10} Minimum Singular Value')
        
        hold off
        
        % % Plot the condition numbers of each subresultant
        figure_name = sprintf('%s : Condition S_{k}',mfilename);
        figure('name',figure_name)
        title('Condition of each subresultant S_{k,k}')
        hold on
        plot(log10(vCondition),'-s','DisplayName','Preprocessed');
        %plot(log10(cond_vec_unproc),'-s','DisplayName','Unprocessed');
        xlabel('k : index of subresultant S_{k}')
        ylabel('log_{10} Condition Number')
        legend(gca,'show')
        hold off
        
        % %
        figure_name = sprintf('%s - Diagonal Norms',mfilename);
        figure('name',figure_name)
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        axis([0.9,min(m,n),-inf,+inf])
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
    case 'n'
    otherwise
        error('error: plot_graphs must be either (y) or (n)')
end