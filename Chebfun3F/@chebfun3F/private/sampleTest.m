function happy = sampleTest(f, cf3Fhandle, tol, enableLog)
% Generate halton sequence
[xeval, yeval, zeval] = halton(30, [-1,1,-1,1,-1,1]);

% Evaluate function
for i = 1:30
    vFun(i) = f(xeval(i),yeval(i),zeval(i));
    vHandle(i) = cf3Fhandle(xeval(i),yeval(i),zeval(i));
end

% Perform tests
if enableLog == 1
    fprintf('Sample Test: error %.2d tol %.2d ', max(abs(vHandle - vFun),[],'all'), 10*tol)
end
fail = (any(max(abs(vHandle - vFun)) > 10*tol));
happy = ~fail;
end


function [x, y, z] = halton(numpts, domain)
% Halton sequences are sequences used to generate points in space, which
% are deterministic and of low discrepancy. They appear to be random for
% many purposes.
%
% From Chebfun2/sampleTest: Adapted from Grady Wright's code 22nd May 2014.

% generate Halton sequences on [0,1]^3:
ndims = 3; % 3D
p = [2 3 5 7 11 13]; % Prime numbers, i.e., bases to be used in 1D Halton sequence generation.
H = zeros(numpts, ndims);
for k = 1:ndims
    N = p(k); v1 = 0; v2 = 0:N-1; lv1 = 1;
    while ( lv1 <= numpts )
        v2 = v2(1:max(2,min(N,ceil((numpts+1)/lv1))))/N;
        [x1,x2] = meshgrid(v2,v1);
        v1 = x1+x2;
        v1 = v1(:);
        lv1 = length(v1);
    end
    H(:,k) = v1(2:numpts+1);
end
% scale [0,1]^3 to the domain of the chebfun3s.
x = H(:,1);
x = (domain(2) - domain(1))*x + domain(1);
y = H(:,2);
y = (domain(4) - domain(3))*y + domain(3);
z = H(:,3);
z = (domain(6) - domain(5))*z + domain(5);
end