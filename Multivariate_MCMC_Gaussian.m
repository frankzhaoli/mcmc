%Multivariate MCMC
%function for analytic value
f=@(x) exp(1).^(-x.^2);

%initialization
sampleSize=5000;
a=-5;
b=5;
r=b-a;
%dimensions>1
dim=4;
mean=0;
sd=1;
start=zeros(1, dim);
%temp to hold samples
temp=zeros(sampleSize, 1);

%target distribution, [chainSize 1] vector
pdf=@(x) mvnpdf(x, zeros(1,dim), eye(dim));
%proposal pdf, [chainSize 1] vector
proppdf=@(x,y) prod(mvnpdf(mean, sd));
%random number generator, [chainSize dim] matrix
proprnd=@(x) unifrnd(a, b, 1, dim); %normrnd(mean, sd, 1, dim);

%gather samples, returns size [sampleSize dim]
sample=mhsample(start, sampleSize, 'logpdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd);

%for each sample
for i=1:sampleSize
    sam=sample(1:i, :);
    temp(i)=prod(sum(f(sam))/i*r);
end
%analytic value
analyticVal=integral(f, a, b)^dim;

%plotting f(x) at each sample
figure;
plot(1:sampleSize, temp);
hold on;
plot([0 sampleSize], [analyticVal analyticVal]);
hold off;
title(['Mean: ', num2str(mean),', SD: ', num2str(sd),', Dimensions: ', num2str(dim)]);
ylabel('f(x)');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(analyticVal))}, 'FontSize', 12, 'TextColor', 'blue');
