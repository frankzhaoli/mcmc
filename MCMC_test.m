<<<<<<< HEAD
%ORIGINAL

%init
H=9.218;
sampleSize=5000;
dim=9;
burn=100;
%target mean/variance
mu=[1, 2, 3, 1, 2, 3, 1, 2, 3];
sigma=[.5, .5, .5, .25, .25, .25, .75, .75, .75];

samples=randn(sampleSize, dim);
samples2=randn(sampleSize, dim);
samples3=randn(sampleSize, dim);
valArray=zeros(sampleSize, 1);
valArray2=zeros(sampleSize, 1);
valArray3=zeros(sampleSize, 1);

%target dis
tpdf=@(x) mvnpdf(x, mu, sigma);

%Metropolis-Hastings
i=1;
while i<sampleSize
    i=i+1;
    %c1
    sample=mvnrnd(samples(i-1, :), sigma);
    a=min([1, tpdf(sample)/tpdf(samples(i-1, :))]);
    
    r=rand();
    if r<a
        samples(i, :)=sample;
    else
        samples(i, :)=samples(i-1, :);
    end
    %c2
    sample=mvnrnd(samples2(i-1, :), sigma);
    a=min([1, tpdf(sample)/tpdf(samples2(i-1, :))]);
    
    r=rand();
    if r<a
        samples2(i, :)=sample;
    else
        samples2(i, :)=samples2(i-1, :);
    end
    %c3
    sample=mvnrnd(samples3(i-1, :), sigma);
    a=min([1, tpdf(sample)/tpdf(samples3(i-1, :))]);
    
    r=rand();
    if r<a
        samples3(i, :)=sample;
    else
        samples3(i, :)=samples3(i-1, :);
    end
end

%Calculate values
for i=1:sampleSize
    sam=samples(1:i, :);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray(i)=real(val);
    
    sam=samples2(1:i, :);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray2(i)=real(val);
    
    sam=samples3(1:i, :);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray3(i)=real(val);
end

%Function value from all samples
val=log(det(2*pi*exp(1)*cov(sam)))/2;
valSum=sum(valArray(burn:sampleSize))/(sampleSize-burn)

total=zeros(sampleSize, 1);
valArray
%Iterating sums
for i=1:sampleSize
    total(i)=sum(valArray(i));
end

%Show histogram/s
%figure;
%hist(samples(:, 3), 30);

%Plot graph
figure;
hold on;
plot(1:sampleSize, total);
%plot(1:sampleSize, valArray2);
%plot(1:sampleSize, valArray3);
plot([0 sampleSize], [H H])
hold off;
ylabel('Entropy');
xlabel('Sample Size');
legend({strcat('Est Val: ', num2str(valSum))}, 'FontSize', 12, 'TextColor', 'blue')

=======

H=9.218;
sampleSize=10000;
i=1;
dim=9;

mu=[1, 2, 3, 1, 2, 3, 1, 2, 3];
sigma=[.5, .5, .5, .25, .25, .25, .75, .75, .75];

samples=zeros(sampleSize, dim);

samples(1, :)=randn(1, dim);

%target
tpdf=@(x) mvnpdf(x, mu, sigma);

%function
f=@(x) mvnpdf(x, mu, sigma);
val=0;

while i<sampleSize
   i=i+1;
   
   sample=mvnrnd(samples(i-1, :), eye(dim));
   a=min([1, tpdf(sample)/tpdf(samples(i-1, :))]);
   
   r=rand();
   if r<a
       samples(i, :)=sample;
   else
       samples(i, :)=samples(i-1, :);
   end
   funcval=f(samples(i, :)).*log(f(samples(i, :)))
   val=val+funcval;
end
val
val/sampleSize


figure;
hist(samples(:, 3), 30);

%{
%plot
figure;
hold on;
plot(1:sampleSize);
plot([0 sampleSize], [H H])
hold off;
%title(strcat('MI: (', num2str(MI), ') vs Sample Size: (', num2str(sampleSize), ')'));
%ylabel('MI');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(H))}, 'FontSize', 12, 'TextColor', 'blue')
%}

%Metropolis-Hastings Algo
function [val]=metHas(sample, prev, mu, sigma)
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
>>>>>>> a8d0114afb57bfa8c448ceb81ab71766409e904a
