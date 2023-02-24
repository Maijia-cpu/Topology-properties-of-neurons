% input: loaded tree using the command: tree=load_tree('tree.mtr');
% output: p : perfection index
%         p_sd: standard devation of perfection index 
%         leaf_N and leaf_num: x and y data points after logarithmic binning
%         error_N: error bar for each data points
tree=load_tree('tree.mtr');
[p,p_sd,leaf_N,leaf_num,error_N]=perfection_index(tree);