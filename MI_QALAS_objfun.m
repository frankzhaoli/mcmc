
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MI, MIimg]=MI_QALAS_objfun(pdv)

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

%% Assign Acquisition Parameters
% Default Parameters
flipAngle=acqparam(1);
TR=acqparam(2);
TE_T2prep=acqparam(3);
Tacq=acqparam(4);
TDpT2=acqparam(5);
TDinv=acqparam(6);
nacq=acqparam(7);
TD=acqparam(8:6+nacq);

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[x,xn,xm,w,wn]=GaussHermiteNDGauss(NumQP,[tisinput(3,1:3),tisinput(5,1:3)],[tisinput(4,1:3),tisinput(6,1:3)]);
lqp=length(xn{1}(:));
parfor qp=1:lqp
    dt=[0,TE_T2prep,Tacq,TDpT2,0,TDinv,Tacq,TD(1),Tacq,TD(2),Tacq,TD(3),Tacq,TD(4)];
    [~,Mmodel_GM(:,qp)]=qalas1p(tisinput(1,1),tisinput(1,1),xn{1}(qp),xn{4}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_WM(:,qp)]=qalas1p(tisinput(1,2),tisinput(1,2),xn{2}(qp),xn{5}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
    [~,Mmodel_CSF(:,qp)]=qalas1p(tisinput(1,3),tisinput(1,3),xn{3}(qp),xn{6}(qp),TR,TE_T2prep,flipAngle,nacq,dt);
end

%% Gauss-Hermite Quadrature MI Approximation
nd=ndims(materialID);
evalstr=sprintf('kspace(jjj%s)=fftshift(fftn(squeeze(imspace(jjj%s))));',repmat(',:',[1,nd]),repmat(',:',[1,nd]));
imspace=0; Ezr=0; Ezi=0; Sigrr=0; Sigii=0; Sigri=0;
materialIDtemp=repmat(permute(materialID,[nd+1,1:nd]),[nacq,ones([1,nd])]);
for iii=1:length(wn(:))
    imspace=zeros(size(materialIDtemp))+(materialIDtemp==1).*Mmodel_GM(:,iii)+(materialIDtemp==2).*Mmodel_WM(:,iii)+(materialIDtemp==3).*Mmodel_CSF(:,iii);
    for jjj=1:size(imspace,1)
        eval(evalstr);
    end
    ksr=real(kspace);
    ksi=imag(kspace);
    
    Ezr=Ezr+ksr*wn(iii);
    Ezi=Ezi+ksi*wn(iii);
    Sigrr=Sigrr+ksr.^2*wn(iii);
    Sigii=Sigii+ksi.^2*wn(iii);
    Sigri=Sigri+ksr.*ksi*wn(iii);
end

N=length(xn);
% std of patient csf = 9.8360; max signal in patient brain = 500;
% max approx signal in synthdata = 0.0584
% std of noise in patient raw data = 17.8574; max signal approx 3000;
signu=3.4762E-4;
% if B1inhomflag==1
%     signu=signu*(1+flipAngle/1.2);
% end
detSigz=(pi^(-N/2)*signu^2 + Sigrr - Ezr.^2).*(pi^(-N/2)*signu^2 + Sigii - Ezi.^2) - (Sigri - Ezr.*Ezi).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MIimg=Hz-Hzmu;

szmi=size(MIimg);
szmi(2:1+ndims(subsmplmask{pdv}))=1;
subsmplmask{pdv}=permute(subsmplmask{pdv},[ndims(subsmplmask{pdv})+1,1:ndims(subsmplmask{pdv})]);
subsmplmask{pdv}=repmat(subsmplmask{pdv},szmi);
MI=sum(MIimg(:).*subsmplmask{pdv}(:));

end