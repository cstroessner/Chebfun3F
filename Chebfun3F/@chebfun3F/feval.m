function out = feval(cf3F, x,y,z)
    u = cf3F.U(x);
    v = cf3F.V(y);
    w = cf3F.W(z);
    
    out = 0;
    r = cf3F.rank();
    for i = 1:r(1)
        for j = 1:r(2)
            for k = 1:r(3)
                out = out + cf3F.C(i,j,k)*u(i)*v(j)*w(k);
            end
        end
    end
end

