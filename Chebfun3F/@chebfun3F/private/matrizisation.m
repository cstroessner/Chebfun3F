function A = matrizisation(T,mode)
s(1) = size(T,1);
s(2) = size(T,2);
s(3) = size(T,3);
if mode == 1
    A = zeros(s(1),s(2)*s(3));
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:s(3)
                A(i,j+s(2)*(k-1)) = T(i,j,k);
            end
        end
    end
elseif mode == 2
    A = zeros(s(2),s(1)*s(3));
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:s(3)
                A(j,i+s(1)*(k-1)) = T(i,j,k);
            end
        end
    end
else
    A = zeros(s(3),s(1)*s(2));
    for i = 1:s(1)
        for j = 1:s(2)
            for k = 1:s(3)
                A(k,i+s(1)*(j-1)) = T(i,j,k);
            end
        end
    end
end
end