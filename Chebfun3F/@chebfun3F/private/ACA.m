%% Adaptive Cross Approximation with full pivoting
function [Ac, At, Ar, rowInd, colInd] = ACA(A, tol, maxIter)

Ac = [];
Ar = [];
At = [];
rowInd = [];
colInd = [];
Aoriginal = A;

for iter = 1:maxIter
    
    [error,I2] = max(abs(A(:)));
    if isempty(error) || error < tol
        return
    end
    
    [I,J] = ind2sub(size(A), I2);
    rowInd = [rowInd, I];
    colInd = [colInd, J];
    
    %A \approx Ac inv(At) Ar'
    Ac = Aoriginal(:,colInd);
    Ar = Aoriginal(rowInd,:)';
    At = Aoriginal(rowInd,colInd);
    
    Anew = zeros(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            Anew(i,j) = A(i,j) - A(i,J)*(1./A(I,J))*A(I,j);
        end
    end
    A = Anew;
end
end