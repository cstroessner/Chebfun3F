clear
testfunctionDataGenerator();
global selectedFunction 
global numfevals 
selectedFunction = 1;

%%
rng(1)
numfevals = 0;
cf = chebfun3(@(x,y,z)testfunction(x,y,z),'vectorize',1);
fprintf('Chebfun3 requires  %10.i function evaluations.\n',numfevals);

%%
rng(1)
numfevals = 0;
cf3F = chebfun3F(@(x,y,z)testfunction(x,y,z), 'enableLog', 0);
fprintf('Chebfun3F requires %10.i function evaluations.\n',numfevals);