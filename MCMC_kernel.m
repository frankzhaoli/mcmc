
%init
load MI_QALAS_objfun_kernel_input.mat;
sampleSize=1000;
dim=9;
dims=1:9;
i=1;
eta=zeros(1, dim);

ksrtemp=zeros(5,151,181);
ksitemp=zeros(5, 151, 181);

%prior mean and sd
pmu=eta_prior(:, 1)';
psigma=eta_prior(:, 2);
%target mean/sd
tmu=eta_target(:, 1);
tsigma=eta_target(:, 2);

samples=zeros(sampleSize, dim);
MItemp=zeros(1, sampleSize);

while i<sampleSize
    %increment
    i=i+1;
    
    for j=1:dim
        %sample from proposal probability distribution
        sample=normrnd(samples(i-1, j), psigma(j));
        samples(i, j)=metHas(sample, samples(i-1, j), tmu(j), tsigma(j));
    end
    eta=abs(samples(i, :));
    [ksr, ksi]=MI_QALAS_objfun_kernel(eta);
    ksrtemp=ksrtemp+ksr;
    ksitemp=ksitemp+ksi;
    
    [MI, MIimg]=calcMI(ksrtemp, ksitemp, i);
    MItemp(i)=MI;
end

imagesc(squeeze(MIimg(1, :, :)));

%plot MI at each sampleSize
figure;
analyticVal=1.1671e+06;
plot(1:sampleSize, MItemp);
hold on;
plot([0 sampleSize], [analyticVal analyticVal])
hold off;
title(['MI vs Sample Size']);
ylabel('MI');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(analyticVal))}, 'FontSize', 12, 'TextColor', 'blue')


%{
for i=1:5
   figure;
   imagesc(squeeze(MIimg(i, :, :)));
end
%}

%{
figure;
for i=1:dim
    subplot(4, 4, i);
    hist(samples(:, i), 20);
    title(strcat('tmu: ', num2str(tmu(i)), ' tsigma: ', num2str(tsigma(i))));
    xlabel(strcat('std: ', num2str(std(samples(:, i))), ' mean: ', num2str(mean(samples(:, i)))));
end
%}

function val=metHas(sample, prev, mu, sigma)
    %target distribution
    f=@(x) normpdf(x, mu, sigma);
    %draw random number
    r=rand();
    %determine acceptance
    alpha=min([1, f(sample)/f(prev)]);
    
    if (r<=alpha)
        %accept proposal
        val=sample;
    else
        %decline proposal
        val=prev;
    end
end

function [MI, MIimg]=calcMI(ksr, ksi, sampleSize)
load MI_QALAS_objfun_kernel_input.mat;
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
szmi(2:1+ndims(subsmplmask{pdv}))=1;
subsmplmask{pdv}=permute(subsmplmask{pdv},[ndims(subsmplmask{pdv})+1,1:ndims(subsmplmask{pdv})]);
subsmplmask{pdv}=repmat(subsmplmask{pdv},szmi);
MI=sum(MIimg(:).*subsmplmask{pdv}(:));
end