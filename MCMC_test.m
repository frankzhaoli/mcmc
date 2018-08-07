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
