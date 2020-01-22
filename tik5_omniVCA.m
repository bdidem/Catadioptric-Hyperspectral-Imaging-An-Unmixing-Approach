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
%%%%% datada ignore edilen siyah ve beyaz pikseller datadan çýkartýlarak
%%%%% endmember estimation yapýlýyor
function [ U ] = tik5_omniVCA(HSI,p,verbose,cozunurluk)
    [m,n,k]=size(HSI.aynaRef);
    M=reshape(HSI.aynaRef,m*n,k)';
    
    [ind_0]=find(all(M==0));
    [ind_1]=find(all(M==1));
    vector=1:size(M,2);
    vector=setdiff(vector,ind_0);
    vector=setdiff(vector,ind_1);
    M_short=M(:,vector);
    
    if nargin==3 
               
        [U(:,:), Paux, snrEstimate ] = hyperVcaModified( M_short, p);
    elseif nargin==4 %OMNI case    
            cozunurluk=cozunurluk(1,vector);

        [U(:,:), Paux, snrEstimate ] = hyperVcaModified( M_short, p, cozunurluk);
    end
        Paux2=vector(Paux);
        % show the pure pixel locations
        locations=(1:(m*n));
        locations=reshape(locations,m,n);
        for i=1:p
            [P(i,1),P(i,2)]=find(locations==Paux2(i));
        end
            

    if strcmp(verbose,'on')
        figure,plot(HSI.wavelength,squeeze(U(:,:))),title('VCA');
               
               
        figure,imshow(HSI.aynaRGB),
        y=squeeze(P(:,2));
        x=squeeze(P(:,1));
        i=(y-1)*m+x;
        hold on,plot(HSI.grid_ind_j(i),HSI.grid_ind_i(i),'r*'),title('Estimated Endmember Locations: VCA ');
        
    end
end