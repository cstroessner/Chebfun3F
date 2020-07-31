rng(1)
Z1 = computeIndicatorConvolution(50000,0.1);
fGP1 = GaussianPeaks(10, 0.5);
fGP2 = GaussianPeaks(3, 0.05);
fGP3 = GaussianPeaks(7, 5);
fGP4 = GaussianPeaks(1, 1);

model = createpde();
geometryFromEdges(model,@squareg);
generateMesh(model,'Hmax',0.25);
%pdemesh(model)
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
fPDE = @(x,y,z) pPDE(x,y,z,model);

save('testfunctionData','Z1','fGP1','fGP2','fGP3','fGP4','fPDE')


function f = GaussianPeaks(numPeaks, avgSigma)
%GAUSSIANPEAKS creates a test function based on several added Gaussian
%density functions
%   numPeaks = number of peaks
%   avgVariance = average covariance of the peaks
f = @(x,y,z) 0;
for iter = 1:numPeaks
    sigma = diag([rand(1)*2*avgSigma,(rand(1)+0.5)*avgSigma,(rand(1)+0.5)*avgSigma]);
    mu = [rand(1)*2-1, rand(1)*2-1, rand(1)*2-1];
    f = @(x,y,z) f(x,y,z) + mvnpdf([x,y,z], mu, sigma);
end
end

function val = indicator1D(x)
val = 0;
if x < 0.5 && x > -0.5
    val = 1;
end
end

function Z1 = computeIndicatorConvolution(n,sig)
t = linspace(-1,1,n);
X = zeros([n,1]);
Y = zeros([n,1]);
g = @(x) normpdf(x,0,sig);
for i = 1:n
    X(i) = indicator1D(t(i));
    Y(i) = g(t(i));
end
Y = Y/sum(Y); %ensure it is a quadrature formula
XF = fft(X,2*n+1);
YF = fft(Y,2*n+1);
Z1 = ifft(YF.*XF);
end

function val = pPDE(p1,p2,p3, model)
f1 = @(location,state)sin(location.y)+cos(location.x)+2;
f2 = @(location,state)cos(location.y)+sin(location.x)+2;
f3 = @(location,state)cos(location.y.^2+location.x.^2)+2;
f = @(location,state) (p1+2)*f1(location,state) + (p2+2)*f2(location,state) + (p3+2)*f3(location,state);

specifyCoefficients(model,'m',0,...
    'd',0,...
    'c',f,...
    'a',0,...
    'f',1);
results = solvepde(model);
%u = results.NodalSolution;
%pdeplot(model,'XYData',u)
val = results.interpolateSolution(0.5,0.5);
end
