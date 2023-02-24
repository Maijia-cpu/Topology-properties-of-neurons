% input: loaded tree using the command: tree=load_tree('tree.mtr');
% output: p : perfection index
%         p_sd: standard devation of perfection index 
%         leaf_N and leaf_num: x and y data points after logarithmic binning
%         error_N: error bar for each data points
function [p,p_sd,leaf_N,leaf_num,error_N]=perfection_index(tree)
    Tchild=child_tree(tree,T_tree(tree));
    BCT_num=typeN_tree(tree);
    Tchild(find(BCT_num==1))=[];
    T0_ind=find(Tchild==0);
    Tchild(T0_ind)=1;% get leaf number for each branch
    clear leafnum leafN 
%% statistics on number of branches and corresponding leaf
    leafN=unique(Tchild);
    for i=1:1:length(leafN)
        leafnum(i)=sum(Tchild==leafN(i));
    end
    
%% logarithmic binning
   clear b bintr leaf_N leaf_num error_N x ydata yfit_RMA
   [b,bintr,leaf_N,leaf_num,error_N]=logarithmic_bin(leafN,leafnum);
   alpha=0.05;
   
   % if the number of data points after binning is smaller than 5,
   % we will fit the power law slope using all data points after the binning.
   % Otherwise, we will not use the first data point for the fitting
   if length(leaf_N)>5 
     x=log(leaf_N(2:end));
     y=log(leaf_num(2:end));
     [b,bintr,bintjm] = gmregress(x,y,alpha);
       if isnan(b(2))
        tbl = table(x', y');
        modelfun = @(b,x) b(1) +b(2)*x;  
        pf=polyfit(x,y,1);
        beta0 = [pf(2) pf(1)]; % Guess values to start with.  Just make your best guess.
        % Now the next line is where the actual model computation is done.
        mdl = fitnlm(tbl, modelfun, beta0);
        table_value=mdl.Coefficients;
        value=table2array(table_value);
        b(2)=value(2,1);
        sd_value=value(2,2);
        else
      end
     n = length(y);
     else
     n = length(leaf_N); 
%      x=log(leaf_N);
%      y=log(leaf_num);
   end   
     t = tinv(1-(alpha/2),n-2);
     p=-b(2)/2;
     if ~exist('sd_value','var')
     sd_value=abs(bintr(2,1)-bintr(2,2))/2/t;
     else
     end
     p_sd=sd_value/2;
  
end