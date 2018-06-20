%Univariate MCMC
clear

%function for analytic value
f=@(x) .5*log(2*pi*exp(1)*(x));

%initialization
sampleSize=10000;
mean=0;
sd=2;
%target distribution
pdf=@(x) normpdf(x, mean, sd+.1);
%proposal pdf
proppdf=@(x, y) normpdf(x, mean, sd);
%random number generator
proprnd=@(x) normrnd(mean, sd);

%gather samples
sample=mhsample(1, sampleSize, 'pdf', pdf, 'proprnd', proprnd, 'proppdf', proppdf);

%initialize temp to hold samples
temp=sample;

%for each sample
for i=1:sampleSize
    %calculate variance of all preceding samples
    sam=sample(1:i, 1);
    temp(i)=f(var(sam));
end
%display expected variance of samples
variance=var(sample);
%analytic value
analyticVal=f(sd^2);

%histogram to show normal distribution of samples
%figure;
%hist(sample, 100);

%plotting f(x) at each sample
figure;
plot(1:sampleSize, temp);
hold on;
plot([0 sampleSize], [analyticVal analyticVal])
hold off;
title(['Mean: ', num2str(mean),', SD: ', num2str(sd),', Proposal: Norm']);
ylabel('f(x)');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(analyticVal))}, 'FontSize', 12, 'TextColor', 'blue')
