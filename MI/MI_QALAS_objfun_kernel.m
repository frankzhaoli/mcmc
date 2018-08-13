
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [ksr,ksi]=MI_QALAS_objfun_kernel(eta_param)

%% Load Input Parameters
load MI_QALAS_objfun_kernel_input.mat;

%% Tissue Properties
% GM/WM/CSF M0/T1/T2 Values
%              GM    WM   CSF  Tumor
T1mean = [1200,  900, 4000, 1200]./1000; % s
T1stdd = [ 100,  100,  200,  150]./1000; % s
T2mean = [ 100,   80, 1000,  110]./1000; % s
T2stdd = [   5,    4,   50,   10]./1000; % s
M0mean = [ 0.9,  0.9,  1.0,  0.9];       % relative intensity
M0stdd = [ .05,  .05,  .05,   .1];       % relative intensity

tisinput=[M0mean;M0stdd;T1mean;T1stdd;T2mean;T2stdd];

%% Default Acquisition Parameters
flipAngle = 4;           % deg
TR = 0.005;              % s
TE_T2prep = 0.100;       % s
Tacq = 0.500;            % s
TDpT2 = 0.4;             % s
TDinv = 0.03;            % s
nacq = 5;
TD = [0.5,0.5,0.5,0.5];          % s

acqparam=[flipAngle,TR,TE_T2prep,Tacq,TDpT2,TDinv,nacq,TD];

%% Generate Quadrature Points for MI Calculation
dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
[~,Mmodel_GM(:)]=qalas1p(eta_param(1),eta_param(1),eta_param(2),eta_param(3),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_WM(:)]=qalas1p(eta_param(4),eta_param(4),eta_param(5),eta_param(6),TR,TE_T2prep,flipAngle,nacq,dt);
[~,Mmodel_CSF(:)]=qalas1p(eta_param(7),eta_param(7),eta_param(8),eta_param(9),TR,TE_T2prep,flipAngle,nacq,dt);

%% Gauss-Hermite Quadrature MI Approximation
nd=ndims(materialID);
evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:)+(materialIDtemp==2).*Mmodel_WM(:)+(materialIDtemp==3).*Mmodel_CSF(:);
for jjj=1:size(imspace,1)
    eval(evalstr);
end
ksr=real(kspace);
ksi=imag(kspace);

end
