figure_name = sprintf('%s : Minimum Singular Values',mfilename)
figure('name',figure_name)
plot(log10(vSingularValues));
hold on

hold off