h = get(0,'children');

nFigures = length(h);


mkdir('Figures/');

for i = 1 : length(h)
    saveas(h(i), ['Figures/figure' num2str(nFigures + 1 - i)], 'fig');
    saveas(h(i), ['Figures/figure' num2str(nFigures + 1 - i)], 'eps');
    saveas(h(i), ['Figures/figure' num2str(nFigures + 1 - i)], 'png');
end