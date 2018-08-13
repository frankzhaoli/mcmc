%MATLAB FUNCTION (MHSAMPLE)

%init
H=9.218;
sampleSize = 5000;
dim = 9;
burn=100;
ent=0;
start = zeros(1, dim);
valArray=zeros(sampleSize, 1);

mu=[1, 2, 3, 1, 2, 3, 1, 2, 3];
sigma=[.5, .5, .5, .25, .25, .25, .75, .75, .75];
delta=.5;
%target distribution, [chainSize 1] vector
pdf=@(x) mvnpdf(x, mu, sigma);
%proposal pdf, [chainSize 1] vector
proppdf=@(x,y) prod(unifpdf(y-x, -delta, delta), 2);
%random number generator, [chainSize dim] matrix
proprnd=@(x) x + rand(1, dim)*2*delta - delta;
%Met-Hast matlab function
sample=mhsample(start, sampleSize, 'pdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd);

%Calculate average value
for i=1:sampleSize
    ent=ent+log(mvnpdf(sample(i, :), mu, sigma))*(-1);
    valArray(i)=ent/i;
end
ent=ent/sampleSize;

%Show histogram
%figure;
%hist(sample(:, 1), 30);

%Plot
figure;
hold on;
plot(1:sampleSize, valArray);
plot([0 sampleSize], [H H]);
hold off;
ylabel('Entropy');
xlabel('Sample Size');
ylim([8 20]);
legend(num2str(ent));