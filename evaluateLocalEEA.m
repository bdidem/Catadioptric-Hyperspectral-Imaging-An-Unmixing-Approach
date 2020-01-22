% temelde tik5_grid_v6 ile ayný
% Usage
%   tik5_grid( exp_no, verbose )
% Inputs
%   exp_no - experiment number to be evaluated
%   verbose - "on" shows the related figures
function error = evaluateLocalEEA(exp_no,verbose)
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
    error=zeros(max(exp_range),1);
    kfAll=zeros(max(exp_range),9);
    for e=exp_range
        e
        [HSI] = tik4_readData_v2(e);
        [m,n,k]=size(HSI.aynaRef);
%             noise_type = 'additive';
%             [w Rn] = estNoise(HSI.short,noise_type,'off');
%             [kfOriginal(e) Ek]=hysimeModified(HSI.short,w,Rn,'off');    % estimate the p
load kfOriginal.mat;
        % gridlerin piksellerini oluþtur
%         grids=create_grid(m,n);%þu an için 3mxn boyutunda, 3=circle sayýsý
        grids=create_grid_v2(m,n);%þu an için (3x3)mxn boyutunda, 9=sub-region sayýsý
%         grids=create_grid_regular(m,n);%þu an için 9mxn boyutunda dikdortgen gridler
        
        kf_start=1;
        for g=1:size(grids,1)
            clear M;
            grid_g=squeeze(grids(g,:,:));
            [grid_ind_i grid_ind_j] = find(grid_g>0);
            
            HSI_tmp = HSI;%HSI_tmp m*n*k boyutunda datayý tutacak. sadece grid içindeki degerler dolu, diðer kýsýmlar siyah
            %HSI_tmp abundance estimation fonksiyonu için kullanýlacak.
            HSI_tmp.aynaRef = bsxfun(@times,HSI.aynaRef,double(grid_g));
             HSI_tmp.aynaRGB(:,:,1)=HSI_tmp.aynaRef(:,:,HSI.rgb(1));
            HSI_tmp.aynaRGB(:,:,2)=HSI_tmp.aynaRef(:,:,HSI.rgb(2));
            HSI_tmp.aynaRGB(:,:,3)=HSI_tmp.aynaRef(:,:,HSI.rgb(3)); 
            
            HSI_short = HSI; %HSI_short grid dýþýndaki boþ alanlarý ignore edecek.
            %HSI_short geometrical approaches için kullanýlacak.
             HSI_short.aynaRGB=HSI_tmp.aynaRGB;
            HSI_short.grid_ind_i = grid_ind_i;
            HSI_short.grid_ind_j = grid_ind_j;
            for i=1:size(grid_ind_i,1)
                M(:,i) = squeeze(HSI_tmp.aynaRef(grid_ind_i(i),grid_ind_j(i),:))';
            end
            f = factor(size(M,2));
            HSI_short.aynaRef = reshape(M',max(f),size(M,2)/max(f),k);
            clear f;
            %% ===============================================================
            % ================ endmember number estimation ===================
            % ================================================================
            [ind_0]=find(all(M==0));
            [ind_1]=find(all(M==1));
            vector=1:size(M,2);
            vector=setdiff(vector,ind_0);
            vector=setdiff(vector,ind_1);
            M_short=M(:,vector);
            noise_type = 'additive';
            if isempty(vector)
                break;
            end
            [w Rn] = estNoise(M_short,noise_type,'off');
            [kfHysime Ek]=hysimeModified(M_short,w,Rn,'off');    % estimate the p
            kfAll(e,g)=kfHysime;
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
    %         tmp=reshape(cozunurluk,m,n);
    %         figure,hist(tmp),title(e);
    %         clear tmp;
            cozunurluk=reshape(cozunurluk,m,n);
          
            for i=1:size(grid_ind_i,1)
                cozunurluk_aux(i) = cozunurluk(grid_ind_i(i),grid_ind_j(i));
            end
            cozunurluk=cozunurluk_aux;
            clear cozunurluk_aux aux minaux maxaux;
            
            [ U(1,:,kf_start:kf_start+kfHysime-1) ] = tik5_omniVCA(HSI_short,kfHysime,verbose,cozunurluk);
%             [ U(1:3,:,kf_start:kf_start+kfHysime-1) ] = tik5_geometricalApproaches_v2(HSI_short,kfHysime,verbose,cozunurluk);

            kf_start=kf_start+kfHysime;
           
        end
        %% endmemberlar kümelenecek
%         for gA=1:3
gA=1;
%             [idx,C] = kmeans(uint8(255*mat2gray(squeeze(U(gA,:,:))')),kfOriginal(e));
            [idx,C] = kmeans(squeeze(U(gA,:,:))',kfOriginal(e));

            nit=10;
            err_min=100;
            err_ind=0;
            errorMapMin=ones(m,n)*100;
            for it=1:nit
    %             it
                while 1
                    for c=1:kfOriginal(e)
                        idx_c=find(idx==c);
                        cluster_i(it,c)=idx_c(randi([1 size(idx_c,1)],1)); 
                    end
                    if it>1
                        if sum(ismember(cluster_i(1:it-1,:),cluster_i(it,:),'rows'))==0
                            break;
                        end
                    else
                        break;
                    end
                end

                U_tmp = squeeze(U(gA,:,cluster_i(it,:)));

                    [m,n,k] = size(HSI.aynaRef);
                    M = reshape(HSI.aynaRef,m*n,k)';

                    [ X ] = hyperNnls(M,U_tmp(:,:));
                    %M-UX
                    B=U_tmp(:,:)*X;
                    % error map oluþturulacak
                           tmpM=reshape(M,m*n*k,1);
                            tmpB=reshape(B,m*n*k,1);
                            a=(tmpM.^2-tmpB.^2);
                            a=reshape(a,k,m*n);
                            errorMap=sqrt(mean(a));
    %                         N=(m*n)-size(find(all(M==0)),2);
    %                         error=sum(errorMap)/N;
                            errorMap=reshape(errorMap,m,n);                    
                for m_i=1:m
                    for n_i=1:n

                        if errorMap(m_i,n_i)<errorMapMin(m_i,n_i)
                            errorMapMin(m_i,n_i)=errorMap(m_i,n_i);
                        end
                    end
                end
                clear U_tmp;
            end
            if strcmp(verbose,'on')
                figure,imshow(errorMapMin,[]),title('VCA Error Map');
            end
            errorMapMin=reshape(errorMapMin,1,m*n);
            N=(m*n)-size(find(all(errorMapMin==0)),2);

            error(e,gA)=sum(errorMapMin)/N;
%         end
        %%
        if strcmp(verbose,'on')
            parentFolder=strcat('D:\Google Drive\TEZ\TÝK5\evaluateLocalEEA\circlar regions_9_v2\exp',num2str(e));
            %newFolder=strcat(parentFolder,'\kf',num2str(kf));
            mkdir(parentFolder);
            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 1:length(FigList)
              FigHandle = FigList(iFig);
    %           FigName   = get(FigHandle, 'Name');
                switch iFig
                    case 19
                        FigName = '1_Omni_Endmembers1.fig';
                    case 18
                        FigName = '2_Omni_EMlocations1.fig';
                    case 17
                        FigName = '3_Omni_Endmembers2.fig';
                    case 16
                        FigName = '4_Omni_EMlocations2.fig';
                    case 15
                        FigName = '5_Omni_Endmembers3.fig';
                    case 14
                        FigName = '6_Omni_EMlocations3.fig';
                    case 13
                        FigName = '7_Omni_Endmembers1.fig';
                    case 12
                        FigName = '8_Omni_EMlocations1.fig';
                    case 11
                        FigName = '9_Omni_Endmembers2.fig';
                    case 10
                        FigName = '10_Omni_EMlocations2.fig';
                    case 9
                        FigName = '11_Omni_Endmembers3.fig';
                    case 8
                        FigName = '12_Omni_EMlocations3.fig';
                    case 7
                        FigName = '13_Omni_Endmembers1.fig';
                    case 6
                        FigName = '14_Omni_EMlocations1.fig';
                    case 5
                        FigName = '15_Omni_Endmembers2.fig';
                    case 4
                        FigName = '16_Omni_EMlocations2.fig';
                    case 3
                        FigName = '17_Omni_Endmembers3.fig';
                    case 2
                        FigName = '18_Omni_EMlocations3.fig';
                    case 1
                        FigName = '19_errorMaps.fig';
                end
              savefig(FigHandle, fullfile(parentFolder, FigName));
            end
%             for iFig = 1:length(FigList)
%               FigHandle = FigList(iFig);
%     %           FigName   = get(FigHandle, 'Name');
%                 switch iFig
%                     case 7
%                         FigName = '1_Omni_Endmembers1.fig';
%                     case 6
%                         FigName = '2_Omni_EMlocations1.fig';
%                     case 5
%                         FigName = '3_Omni_Endmembers2.fig';
%                     case 4
%                         FigName = '4_Omni_EMlocations2.fig';
%                     case 3
%                         FigName = '5_Omni_Endmembers3.fig';
%                     case 2
%                         FigName = '6_Omni_EMlocations3.fig';
%                     case 1
%                         FigName = '7_errorMaps.fig';
%                 end
%               savefig(FigHandle, fullfile(parentFolder, FigName));
%             end
        end
%         save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateLocalEEA\circlar regions_18\exp',num2str(e),'\allEndmembers.mat'),'U');  
       
        close all;
        clear U;
        clear cluster_i;
    end
     save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateLocalEEA\circlar regions_9_v2\error17_2.mat'),'error');  
        save(strcat('D:\Google Drive\TEZ\TÝK5\evaluateLocalEEA\circlar regions_9_v2\kfAll17_2.mat'),'kfAll');  
end