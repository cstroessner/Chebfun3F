function [val] = testfunction(x,y,z)
global selectedFunction
switch selectedFunction
    case 1
        val = 1./(1+25*(x.^2+y.^2+z.^2));
    case 2
        val = exp(x.*y.*z);
    case 3
        val = 0;
    case 4
        val = 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) - 1/3*exp(-(x+1).^2 - y.^2);
    case 5
        val = -cos(50*pi*(x+y+z));
    case 6
        val = exp(-(x-y));
    case 7
        val = 1./cosh(3*(x+y+z)).^2;
    case 8
        val = tanh(5*(x+z)).*exp(y);
    case 9
        val = (x.^2-y.^3+1/8);
    case 10
        val = x.^10*y.^10*z.^10;
    case 11
        val = 1;
    case 12
        val = -cos(50*pi*(x+y+z));
    case 13
        val = x.*y+cos(x.*y);
%     case 14 %TODO out of memory 
%         val = 1/(0.00001+x.^2); 
    case 15
        val = -1000+x.*y.*z.^100;
    case 16
        val = y.*z+cos(z.*y);
    case 17
        t = 5;
        val =  exp(-((x-0.6*cos(pi*t)).^2 + (y-0.6*sin(pi*t)).^2)./0.1);
    case 18
        const = (sqrt(5) + 1) / 2;
        val = 8*(x.^2 - const^4*y.^2) .* (y.^2 - const^4*z.^2) .* ...
            (z.^2 - const^4.*x.^2) .* (x.^4 + y.^4 + z.^4 - 2*x.^2.*y.^2 -...
            2*x.^2.*z.^2 - 2*y.^2.*z.^2) + (3 + 5*const) .* ((x.^2 + y.^2 ...
            + z.^2)).^2 .* ((x.^2 + y.^2 + z.^2 - (2-const))).^2;
    case 19
        val = besselj(0,10*sqrt(x.^2+y.^2+z.^2));
    case 20
        val = ((x-0.4).^2 + y.^2) .* ((x+0.4).^2 + y.^2) - z.^4;
    case 21
        val = 81*(x.^3 + y.^3 + z.^3) - 189*(x.^2.*y + x.^2.*z +...
            y.^2.*x + y.^2.*z + z.^2.*x + z.^2.*y) + 54*(x.*y.*z) + ...
            126*(x.*y + x.*z + y.*z) - 9*(x.^2 + y.^2 + z.^2) - ...
            9*(x + y + z) + 1;
    case 22
        val = abs(((x+1i*y).*exp(1i*z)).^2-1).^2;
    case 23
        val = cos(pi + 5 * (x+y+z));
    case 24
        val = 1./(16 + 5 * (x+y+z));
    case 25
        val = exp(-25*(x.^2+y.^2+z.^2));
    case 26
        val =  1./((0.04+x.^2) .* (0.04+y.^2) .* (0.04+z.^2));
    case 27
        val = x.^2/4000 + y.^2/4000 + z.^2/4000 - ...
            cos(x).*cos(y/sqrt(2)).*cos(z/sqrt(3)) + 1;
    case 28
        val =  z.^2 + 0.5*(x.^2+y.^2-.5).^2;
    case 29
        val =  x.^4 + y.^4 + z.^4 - 0.1*(x.^2 +y.^2 + z.^2) ...
            - 0.5*(x.^2 .* y.^2 + x.^2.*z.^2 + y.^2.*z.^2) + 0.1*x.*y.*z - 1;
    case 30
        val =  cos(2*pi*x).^2 + cos(2*pi*y).^2 + cos(2*pi*z).^2;
    case 31
        val = (sin(pi*(1 + (x - 1)/4))).^2 + ((z - 1)/4).^2 .* ...
            (1+(sin(2*pi*(z - 1)/4)).^2) + (((x - 1)/4)).^2 .* ...
            (1+10*(sin(pi*(1 + (x - 1)/4)+1)).^2) + ...
            (((y - 1)/4)).^2 .* (1+10*(sin(pi*(1 + (y - 1)/4)+1)).^2);
    case 31
        val =  100./(1+100*(x.^2+y.^2+z.^2));
    case 32
        val =  30 + (x.^2-10*cos(2*pi*x)) + ...
            (y.^2-10*cos(2*pi*y)) + (z.^2-10*cos(2*pi*z));
    case 33
        val = 100*(y-x.^2).^2 + (x-1).^2 + 100*(z-y.^2).^2 + (y-1).^2;
    case 34
        val = 1./(1+x.^2+y.^2+z.^2);
    case 35
        val = (cos(2*x+1) + 2*cos(3*x+2) + 3*cos(4*x+3) + ...
            4*cos(5*x+4) + 5*cos(4*x+5)) .* (cos(2*y+1) + 2*cos(3*y+2) + ...
            3*cos(4*y+3) + 4*cos(5*y+4) + 5*cos(4*y+5)) .* ...
            (cos(2*z+1) + 2*cos(3*z+2)+ 3*cos(4*z+3) + 4*cos(5*z+4) + ...
            5*cos(4*z+5));
    case 36
        val =  3*x.^7.*z + y.*z + y.*z.^2 + log(2+y).*z.^3 ...
            - 2*z.^5;
    case 37
        val =  exp(sin(50*x)) + sin(60*exp(y)).*sin(60*z) + ...
            sin(70*sin(x)).*cos(10*z) + sin(sin(80*y)) - sin(10*(x+z)) +...
            (x.^2 + y.^2 + z.^2)/4;
    case 38
        val = log(1+x.^2+y.^2+z.^2);
    case 39
        val =  y.*z+cos(z.*y);
    case 40
        val =  10000./(1+10000*(x.^2+y.^2+z.^2));
    case 41
        val = log(x+y.*z+exp(x.*y.*z)+cos(sin(exp(x.*y.*z))));
    case 42
        val =  cos(2*pi*x).^2 + cos(2*pi*y).^2 + cos(2*pi*z).^2;
    case 43
        val = exp(sin(50*x)) + sin(60.* exp(y)) .* sin(60.*y) + sin(70*sin(x)).*cos(10*z)+sin(sin(80*y))-sin(10*(x+z))+(x.^2+y.^2+z.^2)/4;
    case 44
        val = x.*z + x.^2.*y;
    case 45
        val = sin(1./((0.04+x.^2)).*(0.04+y.^2).*(0.04+z.^2));
%     case 46 %TODO out of memory
%         val = abs((((x+1i*y).*exp(1i*z)).^2-1).^2);
%     case 47
%         val = 1/(cosh(3*(x+y+z)).^2);
    case 48
        val = exp(sin(x+2*y+3*z))+y.*z;
    case 49
        val = exp(1i*pi*x.*y);
    case 50
        val = exp(1i*pi*x.*y.*z);
    case 51
        val =  1./((x+1.1)+(y+1.1)+(z+1.1));
    case 52 
        val = exp(-sqrt((x-1).^2+(y-1).^2+ (z-1).^2));
    case 53 
        load('testfunctionData.mat','Z1');
        convoluted = @(x) interp1(linspace(-2,2,size(Z1,1)),Z1,x, 'pchip');
        val = convoluted(x).*convoluted(y).*convoluted(z);
    case 54
        load('testfunctionData.mat','fGP1');
        val = fGP1(x,y,z);
    case 55
        load('testfunctionData.mat','fGP2');
        val = fGP2(x,y,z);
    case 56
        load('testfunctionData.mat','fGP3');
        val = fGP3(x,y,z);
    case 57 
        load('testfunctionData.mat','fGP4');
        val = fGP4(x,y,z);
    case 58
        val = exp(-sqrt( (x-1).^2+(y-1).^2+(z-1).^2));
    case 'PDE' 
        load('testfunctionData.mat','fPDE');
        val = fPDE(x,y,z);
    otherwise 
        val = 0;    
end
global numfevals 
numfevals = numfevals + numel(val);
end
