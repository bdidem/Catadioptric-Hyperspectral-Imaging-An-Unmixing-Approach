
function error = evaluateSpatialFactor(exp_no,verbose)
    if nargin==0
        exp_range=1:17;
        verbose='on';
    elseif nargin==1
        exp_range=exp_no;
        verbose='on';
    elseif nargin==2
        exp_range=exp_no;
    else
        fprintf('Too many input parameters...');
        quit;
    end
    error=zeros(max(exp_range),9);
    for e=exp_range
        [HSI] = tik4_readData_v2(e);
        [m,n,k]=size(HSI.aynaRef);
        M=reshape(HSI.aynaRef,m*n,k)';
 
        %% ===============================================================
        % ================ endmember number estimation ===================
        % ================================================================
        noise_type = 'additive';
        [w Rn] = estNoise(HSI.short,noise_type,'off');
        [kfHysime Ek]=hysimeModified(HSI.short,w,Rn,'off');    % estimate the p
            %% ===============================================================
            % ================ endmember spectra estimation ==================
            % ================================================================
        %% Omni + geometrical approaches
        f=2;
        a=28.095;%mm
        b=23.4125;%mm
        cozunurluk=cozunurlukHesapla_v5(m,n,f,a,b);

        cozunurluk=reshape(cozunurluk,1,m*n);
%             cozunurluk(find(cozunurluk>0)) = mat2gray(cozunurluk(find(cozunurluk>0)));
        aux=cozunurluk(find(cozunurluk>0));
        minaux=min(aux);
        maxaux=max(aux);
        cozunurluk(find(cozunurluk>0)) = 0.5*((aux-minaux)/(maxaux-minaux))+0.5;
        %%                  
        [U(1,:,:,:)] = tik5_geometricalApproaches_v3(HSI,kfHysime,verbose);
%         [error(e,1:3)] = tik4_abundanceEstimation_v2(HSI,squeeze(U(1,:,:,:)),verbose,'no preprocessing');
        [error(e,1:3)] = tik4_abundanceEstimation_v3(HSI,squeeze(U(1,:,:,:)),verbose,'no preprocessing');

        [HIM_SSPP3D.aynaRef,~] = tik4_sspp(e,HSI.aynaRef,kfHysime);
        HIM_SSPP3D.rgb = HSI.rgb;
        HIM_SSPP3D.wavelength = HSI.wavelength;
        [m2,n2,k2]=size(HIM_SSPP3D.aynaRef);
        M=reshape(HIM_SSPP3D.aynaRef,m2*n2,k2)'; 
        [ind_0]=find(all(M==0));
        [ind_1]=find(all(M==1));
        vector=1:size(M,2);
        vector=setdiff(vector,ind_0);
        vector=setdiff(vector,ind_1);
        HIM_SSPP3D.short=M(:,vector);
        HIM_SSPP3D.shortVector=vector;
        [U(2,:,:,:)] = tik5_geometricalApproaches_v3(HIM_SSPP3D,kfHysime,verbose);
        [error(e,4:6)] = tik4_abundanceEstimation_v3(HSI,squeeze(U(2,:,:,:)),verbose,'SSPP');
        
        [U(3,:,:,:)] = tik5_geometricalApproaches_v3(HSI,kfHysime,verbose,cozunurluk);
        [error(e,7:9)] = tik4_abundanceEstimation_v3(HSI,squeeze(U(3,:,:,:)),verbose,'Omni');
       
        %%
        if strcmp(verbose,'on')
            parentFolder=strcat('D:\Google Drive\TEZ\TÝK5\evaluateSpatialFactor\exp',num2str(e));
            %newFolder=strcat(parentFolder,'\kf',num2str(kf));
            mkdir(parentFolder);
            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 1:length(FigList)
              FigHandle = FigList(iFig);
    %           FigName   = get(FigHandle, 'Name');
                switch iFig
                    case 9
                        FigName = '1_NoPP_Endmembers.fig';
                    case 8
                        FigName = '2_NoPP_EndmemberLocations.fig';
                    case 7
                        FigName = '3_NoPP_ErrorMaps.fig';
                    case 6
                        FigName = '4_SSPP_Endmembers.fig';
                    case 5
                        FigName = '5_SSPP_EndmemberLocations.fig';
                    case 4
                        FigName = '6_SSPP_ErrorMaps.fig';    
                    case 3
                        FigName = '7_Omni_Endmembers.fig';
                    case 2
                        FigName = '8_Omni_EndmemberLocations.fig';
                    case 1
                        FigName = '9_Omni_ErrorMaps.fig';
                end
              savefig(FigHandle, fullfile(parentFolder, FigName));
            end

        end
        save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateSpatialFactor\exp',num2str(e),'\allEndmembers.mat'),'U');  
        save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateSpatialFactor\exp',num2str(e),'\error.mat'),'error');
        close all;
        clear U;
    end
%         save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateSpatialFactor\error.mat'),'error');  

end