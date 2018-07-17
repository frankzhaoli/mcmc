%Multivariate MCMC
%function for analytic value
f=@(x) .5*log(2*pi*exp(1)*(x));

%initialization
sampleSize=5000;
%dimensions>1
dim=10;
mean=0;
sd=1;
start=zeros(1, dim);

%target distribution, [chainSize 1] vector
pdf=@(x) mvnpdf(x, zeros(1,dim), eye(dim));
%proposal pdf, [chainSize 1] vector
proppdf=@(x,y) prod(mvnpdf(mean, sd));
%random number generator, [chainSize dim] matrix
proprnd=@(x) normrnd(mean, sd, 1, dim);

%gather samples, returns size [sampleSize dim]
sample=mhsample(start, sampleSize, 'logpdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd);

%display expected covariance of samples and determinant of covariance
c=cov(sample);
determinant=det(c)
f(determinant)

%initialize temp to hold samples
temp=zeros(sampleSize, 1);

%for each sample
for i=1:sampleSize
    %calculate covariance and determinant of preceding samples
    sam=sample(1:i, :);
    size(sam);
    co=cov(sam);
    d=det(co);
    temp(i)=f(d);
end
%analytic value
sdi=eye(dim)*sd^2;
analyticVal=f(det(sdi));

figure;
hist(sample, 50);

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
