function T = evalTensor(I, J, K, ff)
T = zeros(size(I,2),size(J,2),size(K,2));
for i = 1:size(I,2)
    for j =1:size(J,2)
        for k = 1:size(K,2)
            T(i,j,k) = ff(I(i),J(j),K(k));
        end
    end
end
end