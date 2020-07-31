clear

%creates an error / #evaluations table for Chebfun3 and Chebfun3F

global selectedFunction 
global numfevals 
noTests = 58;

evalTable = zeros([noTests,2]);
errorTable = zeros([noTests,2]);
testPoints = rand([1000,3])*2-1;
x = testPoints(:,1);
y = testPoints(:,2);
z = testPoints(:,3);
fprintf('\n   |      Chebfun3       |      Chebfun3F      |  \n')

for n = 1:noTests
selectedFunction = n;

rng(1);
numfevals = 0;
cf = chebfun3(@(x,y,z)testfunction(x,y,z),'vectorize',1);
evalTable(n,1) = numfevals;
for i = 1:1000
    errorTable(n,1) = max(errorTable(n,1),abs(testfunction(x(i),y(i),z(i))-cf(x(i),y(i),z(i))));
end

rng(1);
numfevals = 0;
cf3F = chebfun3F(@(x,y,z)testfunction(x,y,z));
evalTable(n,2) = numfevals;
for i = 1:1000
    errorTable(n,2) = max(errorTable(n,2),abs(testfunction(x(i),y(i),z(i))-cf3F.feval(x(i),y(i),z(i))));
end

fprintf('%2.i | %9.i  %.2e | %9.i  %.2e | \n',n, evalTable(n,1), errorTable(n,1), ...
    evalTable(n,2), errorTable(n,2))
end
fprintf('\n\n')

%%



