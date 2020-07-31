clear
global selectedFunction 
global numfevals 
numfevals = 0;

%% Figure 1
num = 50;
eps = logspace(0,-5,num);
tau = 10^(-10);
tildetau = tau./eps;
for iter = 1:num
    f = @(x,y,z) 1./(x+y+z+3+eps(iter));
    r(iter) = getRankFromGrid(f, tau);
    n(iter) = getPolynomialDegree(f,tau);
end
clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 400])
loglog(eps,n, 'm')
hold on
loglog(eps,r)
loglog(eps, 28*abs(1./log(1+sqrt(eps))),'r')
loglog(eps, 18+2*abs(log(eps)),'b')
xlabel('$\varepsilon$', 'Interpreter','latex')
leg = legend('max degree for interpolation of fibers', 'max rank of truncated HOSVD','$\mathcal{O}(1/\log(1+\sqrt{\varepsilon}))$', '$\mathcal{O}(|\log(\varepsilon)|)$');
set(leg,'Interpreter','latex');
axis([min(eps),max(eps), 5,1e4])
%print -depsc 'RankDegreeStudy'

%% Figure 2 
clear r
load('Rows.mat'); %contains the 2-mode fibers in Chebfun3
SumV =[];
for i = 1:18
    V{i} = rowsValues{i};
end
for i = 1:18
    SumV = [SumV, V{i}];
    r(i) = rank(SumV);
    s(i) = size(SumV,2);
    size(V{i},2);
end
plot(1:18,r,1:18,s)
clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 370 300])
y = [r; s-r]';
bar(y,'stacked')
plots=get(gca, 'Children');
leg = legend(plots([1,2]),'\# redundant columns of $V_m^{\mathsf{BTD}}$' ,'numerical rank of $V_m^{\mathsf{BTD}}$', 'Location', [0.38 0.79 0.1 0.1]);
set(leg,'Interpreter','latex');
xlabel('$m$','Interpreter','latex')           
%print -depsc 'RedundantFibers'

%% Figure 3
% visualization created directly in latex

%% Figure 4
rng(1)
selectedFunction = 1; % val = 1./(1+25*(x.^2+y.^2+z.^2));
f = @(x,y,z)testfunction(x,y,z); 
cf3F = chebfun3F(f, 'enableLog', 1);
cf3 = chebfun3(f);

%% 
% the numbers are hand copied as they are not stored this detailed
% we ignore the 16 evals in Chebfun3 deciding whether a function is
% vectorizable or not
cf3evals1 = [1000+4913+12167+35937+97336,269658,30];
cf3evals2 = [274625,207668,30];
cf3Fevals1 = [2907+13299+29716+35328,15120,4096+30];
cf3Fevals2 = [50960+56355,9792,4913+30];

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
c = categorical({'Phase 1','Phase 2','Phase 3'});
v1 = [cf3evals1;cf3Fevals1]';
h = bar(c,v1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'BeforeRestart'

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
v2 = [cf3evals2;cf3Fevals2]';
h = bar(c,v2);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'AfterRestart'

%% Figure 5a)
rng(1)
selectedFunction = 58; %val = exp(-sqrt( (x-1).^2+(y-1).^2+(z-1).^2));
f = @(x,y,z)testfunction(x,y,z);
cf3F = chebfun3F(f, 'enableLog', 1);
cf3 = chebfun3(f);

%%
cf3evals = [12870920-30-9544652, 9544652, 30];
cf3Fevals = [2170334-56277-30-44160, 44160, 56277+30];

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
c = categorical({'Phase 1','Phase 2','Phase 3'});
v1 = [cf3evals;cf3Fevals]';
h = bar(c,v1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'F5a'

%% Figure 5b)
rng(1)
selectedFunction = 7; %val = 1./cosh(3*(x+y+z)).^2
f = @(x,y,z)testfunction(x,y,z);
cf3F = chebfun3F(f, 'enableLog', 1);
cf3 = chebfun3(f);

%%
cf3evals = [9354836-30 ,0, 30];
cf3Fevals = [3842596-175560-30, 0, 175560+30];

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
c = categorical({'Phase 1','Phase 2','Phase 3'});
v1 = [cf3evals;cf3Fevals]';
h = bar(c,v1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'F5b'

%% Figure 5c)
rng(1)
selectedFunction = 40; %val =  10000./(1+10000*(x.^2+y.^2+z.^2));
f = @(x,y,z)testfunction(x,y,z);
cf3F = chebfun3F(f, 'enableLog', 1);
cf3 = chebfun3(f);

%% 
%cumulative wrt. restarts
cf3evals = [109269332-30*7-35107576, 35107576, 30*7]; 
cf3Fevals = [1603693-30*5-856254, 856254,30*3+20960];

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
c = categorical({'Phase 1','Phase 2','Phase 3'});
v1 = [cf3evals;cf3Fevals]';
h = bar(c,v1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'F5c'

%% Figure 5d)
testfunctionDataGenerator();
rng(1)
selectedFunction = 'PDE'; 
f = @(x,y,z)testfunction(x,y,z);
cf3F = chebfun3F(f, 'enableLog', 1, 'eps', 1e-9);
cf3 = chebfun3(f, 'output', 1, 'eps', 1e-9);

%% 
cf3evals = [7626-30-1683, 1683, 30]; 
cf3Fevals = [3217-30-125-240,240,30+125];

clf
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 300])
c = categorical({'Phase 1','Phase 2','Phase 3'});
v1 = [cf3evals;cf3Fevals]';
h = bar(c,v1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
legend('Chebfun3','Chebfun3F','Interpreter','latex')
ylabel('$\#$ Function Evaluations','Interpreter','latex')
%print -depsc 'F5d'


%% Figure 6a)
rng(1)
selectedFunction = 1; %val = 1./(1+25*(x.^2+y.^2+z.^2))
f = @(x,y,z)testfunction(x,y,z);
clf
set(gca,'fontsize',10)
ntp = 1000; X = rand([ntp,1]); Y = rand([ntp,1]); Z = rand([ntp,1]);

% sparse grids
no = 0;
num = 6;
err = zeros([num,1]);
addminpoints = logspace(1,7,num);
for i = 1:num 
numfevals = 0;
options = spset('Vectorized', 'on', 'GridType','Chebyshev', 'AbsTol', 10^(-i-2), ...    %Clenshaw-Curtis
    'RelTol',1e-13, 'MaxPoints', 1e8, 'MinPoints', 100+addminpoints(i), 'MinDepth', 4, 'MaxDepth', 10, 'DimensionAdaptive', 'on');
z = spvals(f, 3, [-1 1; -1 1; -1 1], options);
no(i) = numfevals;
ip = spinterp(z, X, Y, Z);
for j = 1: ntp
err(i) = err(i) + (ip(j) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i SG: err %.2d cost %7.i \n', i, err(i), no(i))
end
semilogy(no,err,'m-o')
hold on

% chebfun3
num = 12;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i CF: err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'r-o')

% Chebfun3F
num = 12;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3F(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf.feval(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i mCF err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'b-o')
 
leg = legend('Sparse Grids', 'Chebfun3', 'Chebfun3F');
set(leg,'Interpreter','latex');
xlabel('$\#$ Function Evaluations','Interpreter','latex') 
ylabel('Estimated $\mathcal{L}^2$-Error','Interpreter','latex')
set(figure(1), 'Position', [0 0 370 300])
%print -depsc 'SGa'

%% Figure 6b)
rng(1)
testfunctionDataGenerator();
selectedFunction = 54; %10 Gaussians
f = @(x,y,z)testfunction(x,y,z);
clf
set(gca,'fontsize',10)
ntp = 1000; X = rand([ntp,1]); Y = rand([ntp,1]); Z = rand([ntp,1]);

% sparse grids
no = 0;
num = 8;
err = zeros([num,1]);
addminpoints = logspace(1,3,num);
for i = 1:num 
numfevals = 0;
options = spset('Vectorized', 'on', 'GridType','Chebyshev', 'AbsTol', 10^(-i-2), ...    %Clenshaw-Curtis
    'RelTol',1e-13, 'MaxPoints', 1e8, 'MinPoints', 100+addminpoints(i), 'MinDepth', 4, 'MaxDepth', 10, 'DimensionAdaptive', 'on');
z = spvals(f, 3, [-1 1; -1 1; -1 1], options);
no(i) = numfevals;
ip = spinterp(z, X, Y, Z);
for j = 1: ntp
err(i) = err(i) + (ip(j) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i SG: err %.2d cost %7.i \n', i, err(i), no(i))
end
semilogy(no,err,'m-o')
hold on

% chebfun3
num = 7;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i CF: err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'r-o')

% Chebfun3F
num = 7;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3F(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf.feval(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i mCF err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'b-o')
 
leg = legend('Sparse Grids', 'Chebfun3', 'Chebfun3F');
set(leg,'Interpreter','latex');
xlabel('$\#$ Function Evaluations','Interpreter','latex') 
ylabel('Estimated $\mathcal{L}^2$-Error','Interpreter','latex')
set(figure(1), 'Position', [0 0 370 300])
%print -depsc 'SGb'

%% Figure 6c)
rng(1)
selectedFunction = 41; %val = log(x+y.*z+exp(x.*y.*z)+cos(sin(exp(x.*y.*z))));
f = @(x,y,z)testfunction(x,y,z);
clf
set(gca,'fontsize',10)
ntp = 1000; X = rand([ntp,1]); Y = rand([ntp,1]); Z = rand([ntp,1]);

% sparse grids
no = 0;
num = 7;
err = zeros([num,1]);
addminpoints = logspace(1,7,num); 
for i = 1:num 
numfevals = 0;
options = spset('Vectorized', 'on', 'GridType','Chebyshev', 'AbsTol', 10^(-i-2), ...    %Clenshaw-Curtis
    'RelTol',1e-13, 'MaxPoints', 1e8, 'MinPoints', 100+addminpoints(i), 'MinDepth', 4, 'MaxDepth', 10, 'DimensionAdaptive', 'on');
z = spvals(f, 3, [-1 1; -1 1; -1 1], options);
no(i) = numfevals;
ip = spinterp(z, X, Y, Z);
for j = 1: ntp
err(i) = err(i) + (ip(j) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i SG: err %.2d cost %7.i \n', i, err(i), no(i))
end
semilogy(no,err,'m-o')
hold on

% chebfun3
num = 7;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i CF: err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'r-o')

% Chebfun3F
num = 8;
err = zeros([num,1]);
no = 0;
for i = 1:num
numfevals = 0;
cf = chebfun3F(f, 'eps', 10^(-i-2));
no(i) = numfevals;
for j = 1: ntp
err(i) = err(i) + (cf.feval(X(j),Y(j),Z(j)) - f(X(j),Y(j),Z(j))).^2;
end
err(i) = sqrt(err(i)/ntp);
fprintf('%2.i mCF err %.2d cost %7.i \n', i, err(i), no(i))    
end
semilogy(no,err,'b-o')
 
leg = legend('Sparse Grids', 'Chebfun3', 'Chebfun3F');
set(leg,'Interpreter','latex');
xlabel('$\#$ Function Evaluations','Interpreter','latex') 
ylabel('Estimated $\mathcal{L}^2$-Error','Interpreter','latex')
set(figure(1), 'Position', [0 0 370 300])
%print -depsc 'SGc'


%% required functions
function n = getPolynomialDegree(f,tau)
resolution = 15;
coord = @(i) -1+2*(i-1)/(resolution-1);
n = 0;
for i = 1:resolution
    for j = 1:resolution
        f1 = @(x) f(x,coord(i),coord(j));
        f2 = @(x) f(coord(j),x,coord(j));
        f3 = @(x) f(coord(i),coord(j),x);
        cf1 = chebfun(f1, 'eps', tau);
        cf2 = chebfun(f2, 'eps', tau);
        cf3 = chebfun(f3, 'eps', tau);
        n = max([n, length(cf1), length(cf2), length(cf3)]);
    end
end
end

function r = getRankFromGrid(f, tau)
resolution = 150; 
coord = @(i,n) -cos((i-1).*pi/(resolution-1));
M = zeros([resolution,resolution,resolution]);
M1 = zeros([resolution,resolution*resolution]);
for i = 1:resolution
    for j = 1:resolution
        for k = 1:resolution
            M(i,j,k) = f(coord(i),coord(j),coord(k));
            M1(i,j+resolution*(k-1)) = M(i,j,k);
        end
    end
end
r = max([sum(svd(M1)>tau)]); %f symmetric
end
