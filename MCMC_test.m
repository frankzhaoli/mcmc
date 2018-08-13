%Entropy

%init
H=9.218;
sampleSize=5000;
dim=9;
burn=100;
ent=0;

%target mean/variance
mu=[1, 2, 3, 1, 2, 3, 1, 2, 3];
sigma=[.5, .5, .5, .25, .25, .25, .75, .75, .75];

samples=randn(sampleSize, dim);
valArray=zeros(sampleSize, 1);

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
end

%Calculate average value
for i=1:sampleSize
    ent=ent+log(mvnpdf(samples(i, :), mu, sigma))*(-1);
    valArray(i)=ent/i;
end
ent=ent/sampleSize;

%Show histogram/s
%figure;
%hist(samples(:, 3), 30);

%Plot graph
figure;
hold on;
plot(1:sampleSize, valArray);
plot([0 sampleSize], [H H])
hold off;
ylabel('Entropy');
xlabel('Sample Size');
ylim([8 20]);
legend(num2str(ent));