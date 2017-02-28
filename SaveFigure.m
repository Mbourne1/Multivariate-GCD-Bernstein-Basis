function ret = SaveFigure()


hfigs = get(0, 'children')


index = str2double(input('Enter index of figure you want to save :  ','s'));
if isnan(index) || fix(index) ~= index
  disp('Please enter an integer')
end

filename = input('Enter file name : ','s')
%filename(regexp(filename,'[:]')) = [];
    
fpath = 'Figures/';
    
saveas(hfigs(index), [fpath filename '.fig']) %Matlab .FIG file
saveas(hfigs(index), [fpath filename '.eps']) %Windows Enhanced Meta-File (best for powerpoints)
saveas(hfigs(index), [fpath filename '.png']) %Standard PNG graphics file (best for web)
        
end