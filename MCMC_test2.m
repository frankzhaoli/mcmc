%MATLAB FUNCTION (MHSAMPLE)

%init
H=9.218;
sampleSize = 5000;
dim = 9;
chains=3;
burn=100;
start = zeros(chains, dim);
valArray=zeros(sampleSize, 1);
valArray2=zeros(sampleSize, 1);
valArray3=zeros(sampleSize, 1);

mu=[1, 2, 3, 1, 2, 3, 1, 2, 3];
sigma=[.5, .5, .5, .25, .25, .25, .75, .75, .75];
delta=.5;
%target distribution, [chainSize 1] vector
pdf=@(x) mvnpdf(x, mu, sigma);
%proposal pdf, [chainSize 1] vector
proppdf=@(x,y) prod(unifpdf(y-x, -delta, delta), 2);
%random number generator, [chainSize dim] matrix
proprnd=@(x) x + rand(chains, dim)*2*delta - delta;
%Met-Hast matlab function
sample=mhsample(start, sampleSize, 'pdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd, 'nchain', chains);

for i=1:sampleSize
    sam=sample(1:i, :, 1);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray(i)=real(val);
    
    sam=sample(1:i, :, 2);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray2(i)=real(val);
    
    sam=sample(1:i, :, 3);
    val=log(det(2*pi*exp(1)*cov(sam)))/2;
    valArray3(i)=real(val);
    
end
%Function value from all samples
val=log(det(2*pi*exp(1)*cov(sam)))/2;
valSum=sum(valArray(burn:sampleSize))/(sampleSize-burn)

%Show histogram
%figure;
%hist(sample(:, 1), 30);

%Plot
figure;
hold on;
plot(1:sampleSize, valArray);
plot(1:sampleSize, valArray2);
plot(1:sampleSize, valArray3);
plot([0 sampleSize], [H H]);
hold off;
ylabel('Entropy');
xlabel('Sample Size');
legend({strcat('AV: ', num2str(valSum))}, 'FontSize', 12, 'TextColor', 'blue');
