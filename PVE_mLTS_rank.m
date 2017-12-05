%%%%%%%%  Main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data:(1) Brain mask with the same size as ASL data (to remove background)
%                Here, as an example, matrix size: 46x55x46, it could be
%                any different size
%            (2) Grey matter (GM) mask with the same size as ASL data
%            (3) White matter (WM) mask with the same size as ASL data
%            (4) ASL data

%Output data:(1) PVE-corrected deltaM in GM
%            (2) PVE-corrected deltaM in WM



%Note:(1) To read and write image in analyze format, SPM and Resting-State fMRI Data Analysis Toolkit (REST)should be
%         installed beforehand
%     (2) This program can be applied to a mean CBF map (deltaM in the code) as well as ASL fMRI time series

% Reference: The Matlab code is based on the following paper:
% Please cite: Xiaoyun Liang, Alan Connelly, Fernando Calamante. Improved partial volume correction for single inversion time
%            arterial spin labeling data.Magn Reson Med. 2013 Feb;69(2):531-7. doi: 10.1002/mrm.24279. Epub 2012 Apr 23.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
[mask, header]=rest_ReadNiftiImage('/home/imagetech/DATA/ICA_SWN/AAL/wMNI152lin_T1_1mm_brain.nii');   % Brain mask

[GM_mask, header]=rest_ReadNiftiImage('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/wGM_mask.nii');       % Mask of GM

% [WM_mask, header]=rest_ReadNiftiImage('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/wWM_mask.nii.nii');   % mask of WM
[WM_mask, header]=rest_ReadNiftiImage('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/wWM4.nii');

[Pgm,header]=rest_ReadNiftiImage('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/wGM4.nii');                % Probability map of GM

Pwm=rest_ReadNiftiImage('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/wWM4.nii');                         % Probability map of WM

% Initialization,as an example, a MNI template is used with voxel size:4*4*4, but actually any
% other data can be used.
deltaM_GM=zeros(46,55,46);
deltaM_WM=zeros(46,55,46);
deltaM=zeros(46,55,46,59);
sumGM=zeros(46,55,46);
sumWM=zeros(46,55,46);
meanGM=zeros(46,55,46);
meanWM=zeros(46,55,46);
deltaM_series=zeros(46,55,46,59);   % deltaM is the magnetization difference



%Read ASL fMRI time series
for i=1:1    % i: number of time points in ASL fMRI time series, i=1 for mean deltaM
    if i<10
       filename=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/wrrgrase_im_single00%d.hdr',i);
    else
        if i<100
           filename=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/wrrgrase_im_single0%d.hdr',i);
        else
           filename=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/wrrgrase_im_single%d.hdr',i); 
        end
    end
    [I_ASL,header]=rest_ReadNiftiImage(filename);
    deltaM_series(:,:,:,i)=I_ASL;

end





%%%PVE correction

mag=zeros(46,55,46,59,2);  
mag_LS=zeros(46,55,46,59,2);
kernel_size=zeros(46,55,46);



for u=1:1 %%%   u: number of time points, u=1 for mean deltaM
   if u<10
        filename2=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/ASL_GM_PVC0%d',u);  % Output filename for GM
        filename3=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/ASL_WM_PVC0%d',u);  % Output filename for WM

    else
        filename2=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/ASL_GM_PVC%d',u); % Output filename for GM
        filename3=sprintf('/home/imagetech/DATA/ASL_BOLD_SWN/ASL/WZ/PVC/ASL_WM_PVC%d',u); % Output filename for WM

   end

% Please refer to the MRM papaer
for i=3:44
    for j=3:53
        for k=1:46
            num=0;
            P=zeros(1,2);
            P_rank1=zeros(1,2);
%             
            
            M_rank1=0;
%           
            M=0;
            difference=zeros(5,5);
            P_gm=zeros(5,5);
            P_wm=zeros(5,5);
            M_temp=zeros(5,5);
            if mask(i,j,k)>0.05
                 for m=1:5
                    for n=1:5
                      if mask(i+m-3,j+n-3,k)>0
                        

                           %Sort the voxels according to the difference of the central voxel
                           %from all neighboring voxels in ascending order
                            difference(m,n)=abs(deltaM_series(i,j,k,u)-deltaM_series(i+m-3,j+n-3,k,u));
                            difference_reshape=reshape(difference,25,1);
                            P_gm(m,n)=Pgm(i+m-3,j+n-3,k);
                            P_wm(m,n)=Pwm(i+m-3,j+n-3,k);
                            P_gm_reshape=reshape(P_gm,25,1);
                            P_wm_reshape=reshape(P_wm,25,1);
                            M_temp(m,n)=deltaM_series(i+m-3,j+n-3,k,u);
                            M_temp_reshape=reshape(M_temp,25,1);

                         
                           num=num+1;
                           P=[P;Pgm(i+m-3,j+n-3,k) Pwm(i+m-3,j+n-3,k)];
                           M=[M;deltaM_series(i+m-3,j+n-3,k,u)];
%                        end
                          
                      end
                    end
                 end
                
                 [distance,index]=sort(difference_reshape,'ascend');
                
                  w=0;
                  while(rank(P_rank1)<2&&w<25)
                     w=w+1;
                     P_rank1=[P_rank1;P_gm_reshape(index(w)) P_wm_reshape(index(w))];
                     M_rank1=[M_rank1;M_temp_reshape(index(w))];
                      
                      
                  end
                

            end
            r=rank(P(2:num+1,:));
            kernel_size(i,j,k)=r;
            
            if r<2
                mag(i,j,k,u,:)=[0 0]';
            else
               %%% regaular linear regression
               mag_LS(i,j,k,u,:)=inv((P(2:num+1,:))'*P(2:num+1,:))*(P(2:num+1,:))'*M(2:num+1);
               
               % remove unrealistic values caused by ill conditioned matrix inverse; thresholds for deltaM need to be
               % adjusted accordingly
               if mag_LS(i,j,k,u,1)>80   
                   mag_LS(i,j,k,u,1)=0;
               end
               if mag_LS(i,j,k,u,2)>40
                   mag_LS(i,j,k,u,2)=0;
               end
                   

               % solving the problem with least trimmed squares
               mag(i,j,k,u,:)=mlts_prior(P(2:num+1,:),M(2:num+1),P_rank1,M_rank1,0.4);
                  
               % remove unrealistic values (too high or too low), thresholds for deltaM need to be
               % adjusted accordingly, eg. 40 for WM and 80 for GM
               %for GM
               if mag(i,j,k,u,1)>80 ||mag(i,j,k,u,1)<0
                  mag(i,j,k,u,1)=mag_LS(i,j,k,u,1);
               end
                  
               % for WM   
               if mag(i,j,k,u,2)>40
                  mag(i,j,k,u,2)=mag_LS(i,j,k,u,2);
               end

                  
                  


            end
%             
        end
    end
end

 deltaM1(:,:,:,u)=mag(:,:,:,u,1);   %magnetization difference of GM
 deltaM2(:,:,:,u)=mag(:,:,:,u,2);   %magnetization difference of WM
% 


  sumGM=sumGM+deltaM1(:,:,:,u);
  sumWM=sumWM+deltaM2(:,:,:,u);
  deltaM_GM(:,:,:,u)=deltaM1(:,:,:,u).*GM_mask;
  deltaM_WM(:,:,:,u)=deltaM2(:,:,:,u).*WM_mask;
  
  % Write deltaM of GM
  [flag]=rest_WriteNiftiImage(deltaM_GM(:,:,:,u),header,filename2);
  % Write deltaM of GM
  [flag]=rest_WriteNiftiImage(deltaM_WM(:,:,:,u),header,filename3);
%  
end

