# SpectralPA
spectral PA processing code
% run_sPA_LLS
% CLB V1 Nov 2014 cleaned up linear least squares processing
% parts from Geoff's code
% V3 removed pixel normalization before LLS
% V5 combine hemoglobin into 1 spectra
% V7 Adapted fro embryo imaging session
%%
close all
clear all
clc

%/// user-defined \\\
PAThres = 0;
stepO2 = 0.1;
%load files
load EnergyRatio %energy output from LAZR measured with independent power meter, in a column
load HbSpec %spectra to be fit, first column is oxy, second is deoxy

%initialize matrixes

Ninc = (1/stepO2+1);
Nconc = size(M,2)-1;

pfiles = dir('*sPA.mat');
NFiles = size(pfiles,1); % use when processing multiple spectral PA sessions at a time

%%
for ifile = 1:NFiles
%for ifile = 1:1 %for debugging
    fname = pfiles(ifile).name;
    load(fname)
    fnamesave = pfiles(ifile).name(1:end-8);
    
    %create dimensions
    PA = permute(PA, [1,2,4,3]);
    NWavelength = size(PA,4);
    NFrames = size(PA,3);
    Nj = size(PA,1);
    Nk = size(PA,2);
    
    %prepare file for writing data
    m = matfile([num2str(fnamesave) '_LLS_V7_Conc'], 'Writable', true);
    %m.C_seg = zeros(Nj,Nk,NFrames,Nconc,Ninc);
   %m.S_seg = zeros(Nj,Nk,NFrames,Ninc);
    NsPA = NFrames * Nj * Nk;
    
    %initialize matrixes
    RsPA = zeros(NsPA, NWavelength);
    repEnergyRatio = zeros(NWavelength, NsPA);
    m.C_seg = zeros(NsPA,Nconc,Ninc);
    m.S_seg = zeros(NsPA,Ninc);
    pixel_sPA = zeros(NWavelength, NsPA);

    %scale by power meter energy
    RsPA = reshape(PA, [NsPA, NWavelength]);
    repEnergyRatio = repmat(EnergyRatio,1,NsPA).';
    pixel_sPA = RsPA.*repEnergyRatio;
       
    %calculate model fit for each %sO2
    for NO2 = 1:1:Ninc
        
        %calculate sO2 spectra
        Minc = M(:,1).*NO2/10+M(:,2).*(1-NO2/10);
        
        %solve system of linear equations y=ax+b where y is conc of Hb or O2Hb,
        % a is M (absorbance spectra), x is PA data
        C = (Minc\(pixel_sPA.')).';
        
        %calculate sum of square residual
        repM = permute((repmat(Minc,1,NsPA)), [2,1]);
        repC = repmat(C,1,NWavelength);
        
        S = sum(((repM.*repC)-pixel_sPA),2).^2;
        
        m.C_seg(:,:,NO2) = C;
        
        m.S_seg(:,NO2) = S;
        
    end
    %%
    %find best fit
    [m.Conc, m.I]  = BestFitConc(m.C_seg, m.S_seg);
    m.Nj = Nj;
    m.Nk = Nk;
    m.NFrames = NFrames;
    m.NWavelength = NWavelength;
    m = matfile([num2str(fnamesave) '_LLS_V7_Conc'],'Writable',false);
    clear C RsPA repEnergyRatio pixel_sPA
end
