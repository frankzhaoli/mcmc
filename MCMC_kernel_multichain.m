
%init
load MI_QALAS_objfun_kernel_input.mat;
sampleSize=25000;
dim=9;
dims=1:9;
i=1; j=1;
eta=zeros(1, dim);
eta2=zeros(1, dim);
eta3=zeros(1, dim);

ksrtemp=zeros(5, 151, 181);
ksitemp=zeros(5, 151, 181);

ksrtemp2=zeros(5, 151, 181);
ksitemp2=zeros(5, 151, 181);
ksrtemp2=unifrnd(-50, 50, 5, 151, 181);
ksitemp2=unifrnd(-50, 50, 5, 151, 181);

ksrtemp3=zeros(5, 151, 181);
ksitemp3=zeros(5, 151, 181);
ksrtemp3=unifrnd(-500, 500, 5, 151, 181);
ksitemp3=unifrnd(-500, 500, 5, 151, 181);

%prior mean and sd
pmu=eta_prior(:, 1)';
psigma=eta_prior(:, 2);
%target mean/sd
tmu=eta_target(:, 1);
tsigma=eta_target(:, 2);

samples=zeros(sampleSize, dim);
samples2=zeros(sampleSize, dim);
samples3=zeros(sampleSize, dim);

MItemp=zeros(1, sampleSize);
MItemp2=zeros(1, sampleSize);
MItemp3=zeros(1, sampleSize);
%infeta=zeros(sampleSize, 9);
infksr=zeros(5, 151, 181);
infksi=zeros(5, 151, 181);

%samples(1, dims)=normrnd(5, 1);
%samples2(1, dims)=normrnd(10, 1);
%samples3(1, dims)=normrnd(15, 1);
        
while i<sampleSize
    %increment
    i=i+1;
    
    for j=1:dim
        %sample from proposal probability distribution

        sample=normrnd(samples(i-1, j), psigma(j));
        samples(i, j)=metHas(sample, samples(i-1, j), tmu(j), 1);
        
        sample2=normrnd(samples2(i-1, j), psigma(j));
        samples2(i, j)=metHas(sample2, samples2(i-1, j), tmu(j), 1);
        
        sample3=normrnd(samples3(i-1, j), psigma(j));
        samples3(i, j)=metHas(sample3, samples3(i-1, j), tmu(j), 1);
    end
    eta=abs(samples(i, :));
    eta2=abs(samples2(i, :));
    eta3=abs(samples3(i, :));
    [ksr, ksi]=MI_QALAS_objfun_kernel(eta);
    [ksr2, ksi2]=MI_QALAS_objfun_kernel(eta2);
    [ksr3, ksi3]=MI_QALAS_objfun_kernel(eta3);
    ksrtemp=ksrtemp+ksr;
    ksitemp=ksitemp+ksi;
    
    ksrtemp2=ksrtemp2+ksr2;
    ksitemp2=ksitemp2+ksi2;
    
    ksrtemp3=ksrtemp3+ksr3;
    ksitemp3=ksitemp3+ksi3;
    
    [MI, ~]=calcMI(ksrtemp, ksitemp, i);
    [MI2, ~]=calcMI(ksrtemp2, ksitemp2, i);
    [MI3, ~]=calcMI(ksrtemp3, ksitemp3, i);
    
    MItemp(i)=MI;
    MItemp2(i)=MI2;
    MItemp3(i)=MI3;
    
    %{
    %save infinity vals
    infeta(i, :)=eta;
    if isinf(MI)
        infksr=ksrtemp;
        infksi=ksitemp;
        %break;
    end
    %}
end

[MI, MIimg]=calcMI(ksrtemp, ksitemp, sampleSize);
infksr;
infksi;
MI
%display img
%figure;
%imagesc(squeeze(MIimg(1, :, :)));

%plot MI at each sampleSize
figure;
analyticVal=1.1671e+06;

hold on;
plot(1:sampleSize, MItemp);
plot(1:sampleSize, MItemp2);
plot(1:sampleSize, MItemp3);
plot([0 sampleSize], [analyticVal analyticVal])
hold off;
title(strcat('MI: (', num2str(MI), ') vs Sample Size: (', num2str(sampleSize), ')'));
ylabel('MI');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(analyticVal))}, 'FontSize', 12, 'TextColor', 'blue')

%display histogram of samples
figure;
for i=1:dim
    subplot(4, 4, i);
    hist(samples(:, i), 20);
    title(strcat('tmu: ', num2str(tmu(i)), ' tsigma: ', num2str(tsigma(i))));
    xlabel(strcat('std: ', num2str(std(samples(:, i))), ' mean: ', num2str(mean(samples(:, i)))));
end

%Metropolis-Hastings Algo
function val=metHas(sample, prev, mu, sigma)
    %target distribution
    tpdf=@(x) normpdf(x, mu, sigma);
    %proposal distribution
    ppdf=@(x, y) normpdf(x, y);
    %draw random number
    r=rand();
    %corrective ratio
    c=ppdf(prev, sample)/ppdf(sample, prev);
    %determine acceptance
    alpha=min([1, tpdf(sample)/tpdf(prev)*c]);
    
    if (r<=alpha)
        %accept proposal
        val=sample;
    else
        %decline proposal
        val=prev;
    end
end

%calculate MI/MIimg
function [MI, MIimg]=calcMI(ksr, ksi, sampleSize)
load MI_QALAS_objfun_kernel_input.mat; %#ok<LOAD>
Er=ksr/sampleSize;
Ei=ksi/sampleSize;
Sr=(ksr.^2)/sampleSize;
Si=(ksi.^2)/sampleSize;
Sri=(ksr.*ksi)/sampleSize;

signu=3.4762E-4;
N=6;
detSigz=(pi^(-N/2)*signu^2+Sr-Er.^2).*(pi^(-N/2)*signu^2+Si-Ei.^2)-(Sri-Er.*Ei).^2;
Hz=0.5.*log((2*pi*2.7183)^2.*detSigz);
Hzmu=0.5.*log((2*pi*2.7183)^2.*signu.^4);
MIimg=Hz-Hzmu;

szmi=size(MIimg);
pdv=1;
szmi(2:1+ndims(subsmplmask{pdv}))=1; %#ok<NODEF>
subsmplmask{pdv}=permute(subsmplmask{pdv},[ndims(subsmplmask{pdv})+1,1:ndims(subsmplmask{pdv})]);
subsmplmask{pdv}=repmat(subsmplmask{pdv},szmi);
MI=sum(MIimg(:).*subsmplmask{pdv}(:));
end