classdef chebfun3F

properties
    U
    V
    W
    C
    
    numEvals
    numRestarts
end

methods
    function cf3F = chebfun3F(varargin)
        cf3F = constructor(cf3F, varargin{:});
    end
end

methods (Access = public)
    varargout = rank(cf3F);
    varargout = degree(cf3F);
    varargout = feval(cf3F,x,y,z);
end

end