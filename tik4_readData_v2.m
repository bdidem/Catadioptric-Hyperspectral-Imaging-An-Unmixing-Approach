%Reads data for tik4_compareAllMethods 
%   compareAllMethods
%       Reads the data
%       Extract the mirror region,
%       If the raw data is radiance, converts to reflectance, ignores dark
%       and saturated pixels
%
% Usage
%   [aynaRef] = tik4_readData(exp_no)
% Inputs
%   exp_no - experiment number to be read
% Outputs
%   HSI - A struct contains:
%       .aynaref: Matrix of endmembers (m x n x k) where k is the number of bands.
%       .rgb bands
%       .wavelength
% v2 sadece simulasyon dosyalarýnýn yeri deðiþti
function [HSI] = tik4_readData_v2(exp_no)
    if exp_no<12 %gerçek görüntüler için
        switch exp_no %gerçek görüntüler için
            case 1
                file_path="D:\TEZ çekimler\TÝK1\2008_12_31_23_37_24";%3
            case 2
                file_path="D:\TEZ çekimler\TÝK1\2009_01_01_01_31_51";%4
            case 3
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_14_13_17";%5
            case 4
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_14_20_12";%6
            case 5
                file_path="D:\TEZ çekimler\TÝK4\2018_07_31_12_17_26";%7
            case 6
                file_path="D:\TEZ çekimler\TÝK4\2018_07_31_13_00_53";%8
            case 7
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_12_09_05";
            case 8
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_13_47_25";
            case 9
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_13_49_47";
            case 10
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_14_09_14";
            case 11
                file_path="D:\TEZ çekimler\TÝK2\2017_11_17_14_16_15";                
        end

            if exist(strcat(file_path,'\raw_rd'))==2
                [D,info]=enviread(char(strcat(file_path,'\raw_rd')));
            elseif exist(strcat(file_path,'\raw_0_rd'))==2
                [D,info]=enviread(char(strcat(file_path,'\raw_0_rd')));
            elseif exist(strcat(file_path,'\raw'))==2
                [D,info]=enviread(char(strcat(file_path,'\raw')));
            elseif exist(strcat(file_path,'\raw_0'))==2
                [D,info]=enviread(char(strcat(file_path,'\raw_0')));
            else
                fprintf('Raw data doesnt exist!');
                quit;
            end

        switch exp_no
            case 1 
                ayna=D(139:277,465:663,:);%3
                wrPx=[58,96]; %kameranýn yanýndaki teflondan seçildi
                brPx=[66,190]; %siyah kaðýttan seçildi
                HSI.rgb=[210,95,55];
                ayna=ayna(:,:,1:360);
                HSI.wavelength=[399.16:1.477:929.4030];
            case 2 
                ayna=D(864:1823,1:982,:);%4
                HSI.rgb=[210,95,55];
                wrPx=[252,750]; %beyaz kaðýttan seçildi
                brPx=[394,766]; %siyah kaðýttan seçildi
                ayna=ayna(:,:,1:360);
                HSI.wavelength=[399.16:1.477:929.4030];
           case 3 
                ayna=D(282:619,157:523,:);%5
                HSI.rgb=[421,191,110];
                wrPx=[108,323]; %beyaz kaðýttan seçildi (daha aydýnlýk)
                wrPx=[216,199]; %teflondan seçildi (daha karanlýk)
                brPx=[263,217]; %siyah kapaktan seçildi
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
            case 4 
                ayna=D(96:348,354:720,:);%6
                HSI.rgb=[421,191,110];
                wrPx=[158,255]; % teflondan seçildi
                brPx=[126,300];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];
            case 5 
                ayna=D(105:226,394:590,:);%7
                wrPx=[85,130];
                brPx=[78,163];
                HSI.rgb=[210,95,55];
                ayna=ayna(:,:,1:360);
                HSI.wavelength=[399.16:1.477:929.4030];     
            case 6 
                ayna=D(308:574,365:680,:);%8
                wrPx=[76,274];
                brPx=[248,154];
                HSI.rgb=[210,95,55];
                ayna=ayna(:,:,1:360);
                HSI.wavelength=[399.16:1.477:929.4030];        
            case 7 
                ayna=D(454:714,132:336,:);
                wrPx=[645-454+1,310-132+1];
                brPx=[560-454+1,185-132+1];
                HSI.rgb=[421,191,110];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
            case 8 
                ayna=D(244:388,340:550,:);%8
                wrPx=[281-244+1,484-340+1];
                brPx=[341-244+1,367-340+1];
                HSI.rgb=[421,191,110];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
            case 9 
                ayna=D(210:365,345:553,:);%8
                wrPx=[239-210+1,512-345+1];
                brPx=[313-210+1,371-345+1];
                HSI.rgb=[421,191,110];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
            case 10 
                ayna=D(732:1100,193:557,:);%8
                wrPx=[986-732+1,512-193+1];
                brPx=[1069-732+1,432-193+1];
                HSI.rgb=[421,191,110];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
            case 11 
                ayna=D(200:582,300:668,:);%8
                wrPx=[411-200+1,601-300+1];
                brPx=[522-200+1,369-300+1];
                HSI.rgb=[421,191,110];
                ayna=ayna(:,:,1:720);
                HSI.wavelength=[398.79:0.739:930.1310];      
        end
        [m,n,k]=size(ayna);
        
        %% ayna RGB
        %     aynaRGB(:,:,1)=ayna(:,:,rgb_bands(1));
        %     aynaRGB(:,:,2)=ayna(:,:,rgb_bands(2));
        %     aynaRGB(:,:,3)=ayna(:,:,rgb_bands(3));
        %     figure,imshow(aynaRGB/10),title('Radiance RGB');
        %%
        %% convert radiance to reflectance
            r=5;
            wr=reshape(ayna(wrPx(1)-r:wrPx(1)+r,wrPx(2)-r:wrPx(2)+r,:),((2*r)+1)^2,k);
            wr=max(wr);
            wrM=reshape(repmat(wr,m*n,1),m,n,k);

            br=reshape(ayna(brPx(1)-r:brPx(1)+r,brPx(2)-r:brPx(2)+r,:),((2*r)+1)^2,k);
            br=min(br);
            brM=reshape(repmat(br,m*n,1),m,n,k);

            aynaRef=(ayna-brM)./(wrM-brM);
        %%
        %% ignore dark or saturated pixels
            M=reshape(aynaRef,m*n,k)';
            [~, col]=find(M>1);
            M(:,unique(col))=ones(k,size(unique(col),1));
            [~, col]=find(M<0);
            M(:,unique(col))=zeros(k,size(unique(col),1));
            HSI.aynaRef=reshape(M',m,n,k);
            
            clear col;
        %% apply a mask
            mask = create_mask(m,n);
            HSI.aynaRef = bsxfun(@times,HSI.aynaRef,double(mask));
%         % write aynaRef in envi format
% %             hdr ve info dosyalarýný manual oluþturdum.
%             filepath_write=strcat('D:\Google Drive\TEZ\TÝK4\compare with other methods\exp',num2str(exp_no));
%             load(strcat(filepath_write,'\info',num2str(exp_no),'.mat'));
%             info.samples=n;
%             info.lines=m;
%             info.bands=k;
%             info.wavelength=HSI.wavelength;
%             enviwrite(HSI.aynaRef,info,strcat(filepath_write,'\aynaRef'),strcat(filepath_write,'\aynaRef.hdr'));
        %% aynaRef RGB
%         aynaRGB(:,:,1)=HSI.aynaRef(:,:,HSI.rgb(1));
%         aynaRGB(:,:,2)=HSI.aynaRef(:,:,HSI.rgb(2));
%         aynaRGB(:,:,3)=HSI.aynaRef(:,:,HSI.rgb(3));
%         
%         figure,imshow(aynaRGB);
    elseif exp_no>11 %simulasyon için
        switch exp_no
            case 12
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\noise0.mat';
            case 13
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\noise1.mat';
            case 14
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\noise2.mat';
            case 15
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\Mixed_noise1.mat';
            case 16
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\Mixed_noise2.mat';
            case 17
                file_path='D:\Google Drive\TEZ\TÝK5\simulasyon veriler\tavan eklendi\Mixed_noise0.mat';
        end
        ODI=load(file_path);
        HSI.aynaRef=ODI.ODI;
        [m,n,k]=size(HSI.aynaRef);
        HSI.rgb=[50,20,10];
        HSI.wavelength = [400,409.390000000000 ,418.780000000000 ,428.170000000000 ,437.560000000000 ,446.950000000000 ,456.340000000000 ,465.730000000000 ,475.120000000000 ,484.510000000000 ,493.900000000000 ,503.290000000000 ,512.680000000000 ,522.070000000000 ,531.460000000000 ,540.850000000000 ,550.240000000000 ,559.630000000000 ,569.020000000000 ,578.410000000000 ,587.800000000000 ,597.190000000000 ,606.580000000000 ,615.970000000000 ,625.360000000000 ,634.750000000000 ,644.140000000000 ,653.530000000000 ,662.920000000000 ,672.310000000000 ,681.700000000000 ,691.090000000000 ,700.480000000000 ,709.870000000000 ,719.260000000000 ,728.650000000000 ,738.040000000000 ,747.430000000000 ,756.820000000000 ,766.210000000000 ,775.600000000000 ,784.990000000000 ,794.380000000000 ,803.770000000000 ,813.160000000000 ,822.550000000000 ,831.940000000000 ,841.330000000000 ,850.720000000000 ,860.110000000000 ,869.500000000000 ,878.890000000000 ,888.280000000000 ,897.670000000000 ,907.060000000000 ,916.450000000000 ,925.840000000000 ,935.230000000000 ,944.620000000000 ,954.010000000000 ,963.400000000000 ,972.790000000000 ,982.180000000000 ,991.570000000000 ,1000.96000000000 ,1010.35000000000 ,1019.74000000000 ,1029.13000000000 ,1038.52000000000 ,1047.91000000000 ,1057.30000000000 ,1066.69000000000 ,1076.08000000000 ,1085.47000000000 ,1094.86000000000 ,1104.25000000000 ,1113.64000000000 ,1123.03000000000 ,1132.42000000000 ,1141.81000000000 ,1151.20000000000 ,1160.59000000000 ,1169.98000000000 ,1179.37000000000 ,1188.76000000000 ,1198.15000000000 ,1207.54000000000 ,1216.93000000000 ,1226.32000000000 ,1235.71000000000 ,1245.10000000000 ,1254.49000000000 ,1263.88000000000 ,1273.27000000000 ,1282.66000000000 ,1292.05000000000 ,1301.44000000000 ,1310.83000000000 ,1320.22000000000 ,1329.61000000000 ,1339 ,1348.39000000000 ,1357.78000000000 ,1414.12000000000 ,1423.51000000000 ,1432.90000000000 ,1442.29000000000 ,1451.68000000000 ,1461.07000000000 ,1470.46000000000 ,1479.85000000000 ,1489.24000000000 ,1498.63000000000 ,1508.02000000000 ,1517.41000000000 ,1526.80000000000 ,1536.19000000000 ,1545.58000000000 ,1554.97000000000 ,1564.36000000000 ,1573.75000000000 ,1583.14000000000 ,1592.53000000000 ,1601.92000000000 ,1611.31000000000 ,1620.70000000000 ,1630.09000000000 ,1639.48000000000 ,1648.87000000000 ,1658.26000000000 ,1667.65000000000 ,1677.04000000000 ,1686.43000000000 ,1695.82000000000 ,1705.21000000000 ,1714.60000000000 ,1723.99000000000 ,1733.38000000000 ,1742.77000000000 ,1752.16000000000 ,1761.55000000000 ,1770.94000000000 ,1780.33000000000 ,1789.72000000000 ,1930.57000000000 ,1939.96000000000 ,1949.35000000000 ,1958.74000000000 ,1968.13000000000 ,1977.52000000000 ,1986.91000000000 ,1996.30000000000 ,2005.69000000000 ,2015.08000000000 ,2024.47000000000 ,2033.86000000000 ,2043.25000000000 ,2052.64000000000 ,2062.03000000000 ,2071.42000000000 ,2080.81000000000 ,2090.20000000000 ,2099.59000000000 ,2108.98000000000 ,2118.37000000000 ,2127.76000000000 ,2137.15000000000 ,2146.54000000000 ,2155.93000000000 ,2165.32000000000 ,2174.71000000000 ,2184.10000000000 ,2193.49000000000 ,2202.88000000000 ,2212.27000000000 ,2221.66000000000 ,2231.05000000000 ,2240.44000000000 ,2249.83000000000 ,2259.22000000000 ,2268.61000000000 ,2278 ,2287.39000000000 ,2296.78000000000 ,2306.17000000000 ,2315.56000000000 ,2324.95000000000 ,2334.34000000000 ,2343.73000000000 ,2353.12000000000 ,2362.51000000000 ,2371.90000000000 ,2381.29000000000 ,2390.68000000000 ,2400.07000000000 ,2409.46000000000 ,2418.85000000000 ,2428.24000000000 ,2437.63000000000 ,2447.02000000000];
    end
    %%
    r2=n*25/1000; %ortadaki maskelenecek alaný data boyutuna göre ayarlar
    HSI.aynaRef((m/2)-r2:(m/2)+r2,(n/2)-r2:(n/2)+r2,:)=0;

    %% short (hysime ve EEA için 0 ve 1 lerin olmadýðý short data.)
   M=reshape(HSI.aynaRef,m*n,k)'; 
   [ind_0]=find(all(M==0));
   [ind_1]=find(all(M==1));
   vector=1:size(M,2);
   vector=setdiff(vector,ind_0);
   vector=setdiff(vector,ind_1);
   HSI.short=M(:,vector);
   HSI.shortVector=vector;
end