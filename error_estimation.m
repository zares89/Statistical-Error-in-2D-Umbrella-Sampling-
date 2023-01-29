clear all;close all;
% 2D Extension of the error estimation method for 1D umbrella sampling simulations developed by  Kästner and Thiel, J. Chem. Phys. 124, 234106 (2006)
% Equation numbers refer to equations in the paper by Kästner and Thiel
% (2006). 
load rep.mat % Load the free energy landscape produced by the Nudged_Elastic_Band code 
rep_mat_all=[-0.01519+7.55+0.096 2.326;rep_mat_all];
rep_mat_all=rep_mat_all(1:7,:);
bins_CV=rep_mat_all;
zdist=7.4:0.1:15.3;
coordall1=1:0.1:4.0;
coordall2=3:0.1:6.0;
coordall3=2:0.1:3.0;
coordcell{1}=zdist;coordcell{2}=coordall1;coordcell{3}=coordall2;coordcell{4}=coordall3;

% Find the nearest points in the actual simulations to the points found in
% the minimum energy path code (Nudged_Elastic_Band)
nearest_id=zeros(size(rep_mat_all,1),2);
for i=1:size(rep_mat_all,1)
    ds=[];
    for j=1:length(coordcell{1})
        ds(j)=abs(rep_mat_all(i,1)-coordcell{1}(j));
    end
    [minval,minid1]=min(ds);
    nearest_id(i,1)=minid1;
    if minid1 <=12
        ds=[];
        for j=1:length(coordcell{2})
            ds(j)=abs(rep_mat_all(i,2)-coordcell{2}(j));
            
        end
        [minval,minid]=min(ds);
        nearest_id(i,2)=minid;
    elseif minid1 >12 && minid1 <= 45
        if rep_mat_all(i,2) > 3
            ds=[];
            for j=1:length(coordcell{3})
                ds(j)=abs(rep_mat_all(i,2)-coordcell{3}(j));
                
            end
            [minval,minid]=min(ds);
            nearest_id(i,2)=minid;
        else
            ds=[];
            for j=1:length(coordcell{4})
                ds(j)=abs(rep_mat_all(i,2)-coordcell{4}(j));
                
            end
            [minval,minid]=min(ds);
            nearest_id(i,2)=minid+31;
        end
    else
        ds=[];
        for j=1:length(coordcell{3})
            ds(j)=abs(rep_mat_all(i,2)-coordcell{3}(j));
            
        end
        [minval,minid]=min(ds);
        nearest_id(i,2)=minid;
    end
end
list_CV=nearest_id;
k_b=0.001985875; % kcal/(mol.K)
T=300;
beta=1/(k_b*T);

mydata1={};
for i=1:size(list_CV,1)

    txtfile=strcat(['hist-all-' num2str(list_CV(i,1)) '-' num2str(list_CV(i,2)) '.txt']);
    mydata1{i}=importdata(txtfile);
    mydata1{i}=mydata1{i}(:,2:3);
    z_top=20;
    z_bot=0;
    del_z=0.1;
    coord_top=7;
    del_coord=0.1;
    coord_bot=0;
    z_length=round(z_top/del_z+1);
    coord_length=round(coord_top/del_coord+1);
    num_cvs=zeros(z_length,coord_length);
    %num_coord=zeros(coord_length,1);
    zprev=0;
    
    IDX = uint32(1:size( mydata1{i},1)); 
    for j=1:z_length
        cprev=0;
        for k=1:coord_length
            %ow_inz_frame=find(ow(:,4)<=zprev+del_z && ow(:,4)>=zprev);
            ind=IDX(mydata1{i}(:,1) >= zprev & mydata1{i}(:,1) < zprev+del_z & mydata1{i}(:,2) >= cprev & mydata1{i}(:,2) < cprev+del_coord );
            
            num_cvs(j,k)=length(ind);
            
          %  zprev=zprev+del_z;
            cprev=cprev+del_coord;
        end
        zprev=zprev+del_z;
    end
    
  %  sum_2d=sum(num_cvs);
  %  sum_2d=sum(sum_2d);
  %  prob_cvs{i}=num_cvs./sum_2d;
    numf_window(i)=size(mydata1{i},1);
   % pdfbiased(i)=pdf(mydata1{i}(:,1),x);
end

z_ax=z_bot:del_z:z_top;
coord_ax=coord_bot:del_coord:coord_top;
 % Calculate the variance of the unbiased mean force at window i
 
meanksiw={};
covarksiw={};
sum=zeros(2,2);
for i=1:size(list_CV,1)
    meanksiw{i}=mean(mydata1{i});
    covarksiw{i}=cov(mydata1{i},1);
    sum=sum+covarksiw{i};
end
avg_cvar_mat=1/size(list_CV,1)*sum;
 % Calculate the variance of covariance matrix and average values of
 % (ksi,eta) in 100 segments
 
 nseg=100;
 
 for i=1:size(list_CV,1)
     meanksi_seg={};
     covarksi_seg={};
     invcovara_allseg=[];
     invcovarb_allseg=[];
     invcovarc_allseg=[];
     invcovard_allseg=[];
     meanksi_allseg=[];
     segnum=round(1/nseg*(size(mydata1{i},1)));
     for j=1:nseg
         meanksi_seg{j}=mean(mydata1{i}((j-1)*segnum+1:j*segnum,:));
         covarksi_seg{j}=cov(mydata1{i}((j-1)*segnum+1:j*segnum,:),1);
         invcovarksi_seg=covarksi_seg{j}^(-1);
         meanksi_allseg=[meanksi_allseg;meanksi_seg{j}];
         invcovara_allseg=[invcovara_allseg;invcovarksi_seg(1,1)];
         invcovarb_allseg=[invcovarb_allseg;invcovarksi_seg(1,2)];
         invcovarc_allseg=[invcovarc_allseg;invcovarksi_seg(2,1)];
         invcovard_allseg=[invcovard_allseg;invcovarksi_seg(2,2)];
     end
     mean_mean_ksi(i)= mean(meanksi_allseg(:,1));
     mean_mean_eta(i)= mean(meanksi_allseg(:,2));
     var_mean_ksi(i)= var(meanksi_allseg(:,1));
     var_mean_eta(i)= var(meanksi_allseg(:,2));
     vara(i)=var( invcovara_allseg);
     varb(i)=var( invcovarb_allseg);
     varc(i)=var( invcovarc_allseg);
     vard(i)=var( invcovard_allseg);
 end
 
% zeta_avg=mean(meanzetaw);
% eta_avg=mean(meanetaw);
% zeta_var=var(meanzetaw);
% vzeta_var=var(varzetaw);
% azeta_var=mean(varzetaw);
% eta_var=var(meanetaw);
% veta_var=var(varetaw);
% aeta_var=mean(varetaw);
% bins_CV=rep_mat_all;
P_cvs=[];
sum=zeros(size(bins_CV,1),1);
for i=1:size(bins_CV,1)
    
    for j=1:size(list_CV,1)
        idx_CV1=find(abs(bins_CV(i,1)-z_ax(:))<0.05);
        idx_CV2=find(abs(bins_CV(i,2)-coord_ax(:))<0.05);
        P_cvs(i,j)=(2*pi)^(-1)*det(covarksiw{j})^(-1/2)*...
            exp(-1/2*(bins_CV(i,:)-meanksiw{j})*covarksiw{j}^(-1)*(bins_CV(i,:)-meanksiw{j})');
       % (2*pi)^(-1)*det(covarksiw{j}*(size(mydata1{j},1)-1))^(-1/2)*...
        %    exp(-1/2*(bins_CV(i,:)-meanksiw{j})*(covarksiw{j}*(size(mydata1{j},1)-1))^(-1)*(bins_CV(i,:)-meanksiw{j})');
    sum(i)=sum(i)+P_cvs(i,j)*numf_window(j);
    end
   % sum(i,j)=sum(i,j)+P_cvs(i,j)*numf_window(j);
end
 p_cvs=[];
 for i=1:size(bins_CV,1) 
     for j=1:size(list_CV,1)
         p_cvs(i,j)=P_cvs(i,j)*numf_window(j)/sum(i);
     end
 end
  % Calculate var(dA/dksi) equation 5
  varforce={};
 for j=1:size(list_CV,1)
     for i=1:size(bins_CV,1)
         varforce{j}(i,1:2)=1/beta^2*[(bins_CV(i,1)-meanksiw{j}(1))^2*vara(j)+covarksiw{j}(1,1)^2*var_mean_ksi(j)+...
             (bins_CV(i,2)-meanksiw{j}(2))^2*varb(j)+covarksiw{j}(1,2)^2*var_mean_eta(j);...
             (bins_CV(i,1)-meanksiw{j}(1))^2*varc(j)+covarksiw{j}(2,1)^2*var_mean_ksi(j)+...
             (bins_CV(i,2)-meanksiw{j}(2))^2*vard(j)+covarksiw{j}(2,2)^2*var_mean_eta(j)];
     end
     
 end
 % Calculate var(dA/dksi) equation 9
  weighted_sumvarf_bin=zeros(size(bins_CV,1),2);
 for i=1:size(bins_CV,1)
     %weighted_sumvarf_bin(i)=0;
     for j=1:size(list_CV,1)
         weighted_sumvarf_bin(i,:)= weighted_sumvarf_bin(i,:)+p_cvs(i,j)^2*varforce{j}(i,:);
         
     end
 end
% store mean of the variance of ksi and eta
varksi=zeros(size(list_CV,1),2);
for i=1:size(list_CV,1)
    varksi(i,1)=covarksiw{i}(1,1);
    varksi(i,2)=covarksiw{i}(2,2);
end
meanvarksi(1)=mean(varksi(:,1));
meanvarksi(2)=mean(varksi(:,2));
% Calculate var(delta(A)) equation 15
 varfree=mean(weighted_sumvarf_bin(:,1))*(bins_CV(end,1)-bins_CV(1,1))*sqrt(meanvarksi(1)*sqrt(2*pi)...
     -2*meanvarksi(1))+mean(weighted_sumvarf_bin(:,2))*(bins_CV(end,2)-bins_CV(1,2))*sqrt(meanvarksi(2)*sqrt(2*pi)...
     -2*meanvarksi(2));
% varfree=mean(weighted_sumvarf_bin)*(norm(bins_CV(end,:)-bins_CV(1,:))*sqrt(det(avg_cvar_mat))*sqrt(2*pi)...
%     -2*det(avg_cvar_mat));

1.96*sqrt(varfree)
