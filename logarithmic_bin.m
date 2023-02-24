function [b,bintr,leaf_N,leaf_num,error_N]=logarithmic_bin(leafN,leafnum)
% logarithic binning
% input: raw data of leaf number n and number of branches that support n leaves
% output: b regression coefficients of RMA fitting
%         bint given confidence intervals for b
%         leaf_N and leaf_num: x and y data points after logarithmic binning
%         error_N: error bar for each data points
%% logarithmic binning
    nmax=max(leafN);
    clear leaf_N leaf_num error_N leaf_sum
  % keep the first data point
    leaf_num(1)=leafnum(1);
    leaf_N(1)=leafN(1);
    error_N(1)=sqrt(sum(leafnum(1))+2);
    leaf_sum(1)=sum(leafnum(1));
% start logarithmic binning after the first data point
    i=2;
    num_left=1.5;
    num_right=2.5;
while num_right<nmax+1
    
    ind_N=find(leafN>=num_left & leafN<num_right);
    leaf_num(i)=sum(leafnum(ind_N))/(num_right-num_left);
    leaf_N(i)=(num_right+num_left)/2;
    error_N(i)=sqrt(sum(leafnum(ind_N))+2)/(num_right-num_left);
    leaf_sum(i)=sum(leafnum(ind_N));
       
    i=i+1;
    num_left=num_right;
    num_right=num_right+0.5*2^(i-1);
end
     %% fit with RMA
    clear x y b bintr bintjm
    ind_l=find(leaf_sum>=9);% threshold, when the number within a certain bin near the tail is smaller than this threshold, that bin is neglected 
    ind_c = circshift(ind_l,1);
    delta_c=ind_l-ind_c;
    ind_delta=find(delta_c~=1);
 
    if length(ind_delta)>1 && ~isempty(ind_delta)
        x=log(leaf_N(1:ind_l(end)));
        y=log(leaf_num(1:ind_l(end)));
        leaf_N=leaf_N(1:ind_l(end));
        leaf_num=leaf_num(1:ind_l(end));
        error_N=error_N(1:ind_l(end));
    else
        x=log(leaf_N(ind_l));
        y=log(leaf_num(ind_l));
        leaf_N=leaf_N(ind_l);
        leaf_num=leaf_num(ind_l);
        error_N=error_N(ind_l);
    end
    alpha=0.05;
    [b,bintr,bintjm] = gmregress(x,y,alpha);
    

end