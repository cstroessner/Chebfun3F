function happy = sampleTest(f, cf3Fhandle, tol, enableLog)


% Evaluate function
[x, y, z] = halton(30);
for i = 1:30
    vFun(i) = f(x(i),y(i),z(i));
    vHandle(i) = cf3Fhandle(x(i),y(i),z(i));
end

% Perform tests
if enableLog == 1
    fprintf('Sample Test: error %.2d tol %.2d ', max(abs(vHandle - vFun),[],'all'), 10*tol)
end
fail = (any(max(abs(vHandle - vFun)) > 10*tol));
happy = ~fail;
end


function [x, y, z] = halton(numpts)
p = haltonset(3,'Skip',1);
H = p(1:30,:);
x = H(:,1);
x = 2*x -1;
y = H(:,2);
y = 2*y -1;
z = H(:,3);
z = 2*z -1;
end