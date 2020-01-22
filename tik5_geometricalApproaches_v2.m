%geometricalApproaches estimates endmembers by using geometrical approaches
%(NFINDR, PPI, VCA)
%
% Usage
%   [ U ] = tik4_geometricalApproaches( aynaRef,p,verbose )
% Inputs
%   aynaRef - HSI data as 3D matrix (m x n x k).
%   p - Number of endmembers to find.
%   verbose - "on" shows the endmember locations and spectra
%   cozunurluk - spatial resolution map for the OMNI experiments
% Outputs
%   U - Matrix of endmembers (3 x k x p) for NFINDR, PPI, VCA respectively.

function [ U ] = tik5_geometricalApproaches_v2(HSI,p,verbose,cozunurluk)
    [m,n,k]=size(HSI.aynaRef);
    M=reshape(HSI.aynaRef,m*n,k)';
    if nargin==3 
        [U(1,:,:), Paux] = NFINDRModified( HSI.aynaRef, p);   
        P(1,:,1)=Paux(:,1);
        P(1,:,2)=Paux(:,2);
        
        [U(2,:,:), Paux] = hyperPpiModified(M, p, k);
        % show the pure pixel locations
        locations=(1:(m*n));
        locations=reshape(locations,m,n);
        for i=1:p
            [P(2,i,1),P(2,i,2)]=find(locations==Paux(i));
        end
        
        [U(3,:,:), Paux, snrEstimate ] = hyperVcaModified( M, p);
        % show the pure pixel locations
        locations=(1:(m*n));
        locations=reshape(locations,m,n);
        for i=1:p
            [P(3,i,1),P(3,i,2)]=find(locations==Paux(i));
        end

    elseif nargin==4 %OMNI case
        [U(1,:,:), Paux ] = NFINDRModified( HSI.aynaRef, p, cozunurluk); 
        P(1,:,1)=Paux(:,1);
        P(1,:,2)=Paux(:,2);
        
        [U(2,:,:), Paux] = hyperPpiModified(M, p, k,cozunurluk);
        % show the pure pixel locations
        locations=(1:(m*n));
        locations=reshape(locations,m,n);
        for i=1:p
            [P(2,i,1),P(2,i,2)]=find(locations==Paux(i));
        end
        
        [U(3,:,:), Paux, snrEstimate ] = hyperVcaModified( M, p, cozunurluk);
        % show the pure pixel locations
        locations=(1:(m*n));
        locations=reshape(locations,m,n);
        for i=1:p
            [P(3,i,1),P(3,i,2)]=find(locations==Paux(i));
        end
        
    end
    
    
    
    
    if strcmp(verbose,'on')
        figure,subplot(1,3,1),plot(HSI.wavelength,squeeze(U(1,:,:))),title('NFINDR');
               subplot(1,3,2),plot(HSI.wavelength,squeeze(U(2,:,:))),title('PPI');
               subplot(1,3,3),plot(HSI.wavelength,squeeze(U(3,:,:))),title('VCA');
               
               
        figure,subplot(3,1,1),imshow(HSI.aynaRGB),
        y=squeeze(P(1,:,2));
        x=squeeze(P(1,:,1));
        i=(y-1)*m+x;
        
        hold on,plot(HSI.grid_ind_j(i),HSI.grid_ind_i(i),'r*'),title('Estimated Endmember Locations: NFINDR ');
       
        subplot(3,1,2),imshow(HSI.aynaRGB),
        y=squeeze(P(2,:,2));
        x=squeeze(P(2,:,1));
        i=(y-1)*m+x;
        hold on,plot(HSI.grid_ind_j(i),HSI.grid_ind_i(i),'r*'),title('Estimated Endmember Locations: PPI ');
        
        subplot(3,1,3),imshow(HSI.aynaRGB),
        y=squeeze(P(3,:,2));
        x=squeeze(P(3,:,1));
        i=(y-1)*m+x;
        hold on,plot(HSI.grid_ind_j(i),HSI.grid_ind_i(i),'r*'),title('Estimated Endmember Locations: VCA ');
        
    end
end