function [ff xx] = TheoryFormFunction(xIn,c)
%
% function [ff xx] = TheoryFormFunction(x)
%
% Estimation of the form function and cross-section area for a single
% sphere based on: 
% 
% Gaunard G.C & Uberall, H., 1983. J.Acoust. Soc. Am. 73(1),pp 1-12.
% (see equation 17 in page 3 of paper). 
%
%
%  INPUTS 
%       xIn  = ka where k=wavenumber and a sphere radius
%       c  = speed of sound in water
%  OUTPUTS
%       ff  = Acoustic backscatterance form function (monostatic) 
%      xx = cross-sectional area of the spheres.
%
%  NOTE: The results depend on the property values used for the spheres and
%        the density in particular. Default values are for glass spheres
%        with density=2586kg/m3, shear velocity=3545m/s, compressional
%        velocity=5550m/s. See code for values of Tugsten Carbide and Brass
%        Spheres. Water density is used 1000kg/m and speed of sound is set 
%        to 1476m/s. 
%        These values are those reported in  
%        Thorne P.D. & S.C. Cambell, 1992. Backscattering by suspension of Spheres.
%        J. Acoust. Soc. Am. 92(2):978-986.
%  
% George Voulgaris, 
% CPSD Lab, July, 2007
% gvoulgaris@geol.sc.edu
% Version History
% 1.0 Provided by George Volgaris
% 1.1 Modified by Aquatec to suit the calling function and renamed
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Environmental Properties
rho_w = 1000;              % Density of water at 18C
c_1   = c;                 % Default speed of sound in medium (water, m/s)=1480 
%
% Material/Sphere Properties
% 
% rho_s = 13800;              % Density of Tugsten Carbide
rho_s =  2586;                % Density of Glass Sphere 
% rho_s =  8360;              % Density of Brass Sphere 
rhor=rho_w/rho_s;           % Ratio of densities
%
% Shear wave velocity
% 
%c_s   = 4185;  % Default shear wave velocity for Tugsten Carbid (m/s)
c_s   = 3545;  % Default shear wave velocity for Glass Sphere (m/s)
%c_s   = 2100;  % Default shear wave velocity for Brass (m/s)
%
% Compressional wave velocity
%
%c_d   = 6860;  % Default compressional wave velocity for Tugsten Carbid (m/s)
c_d   = 5550;  % Default compressional wave velocity for Glass Sphere (m/s)
%c_d   = 4372;  % Default compressional wave velocity for Brass (m/s)
%
ff = zeros(size(xIn));
xx = zeros(size(xIn));

for index = 1:size(xIn,2)

    ka = xIn(:,index);    % x=ka for mean radius
    aa = 0.032;    % radius of sphere to use for calculations
    %
    for I=1:length(ka)            % For each x=ka 
      x     = ka(I);              % x=ka value for a(i)
      k_1   = x/aa;               % wavenumber=omega/c_1;
      k_s2  = (c_1/c_s)*k_1;      % omega/c_s;
      k_d2  = (c_1/c_d)*k_1;      % omega/c_d;
      xs2   = (k_s2/k_1)*x;
      xd2   = (k_d2/k_1)*x; 
      %
      % Initialize summation loop
      %
      the_sum=0;                  % for form function
      the_sum1=0;                 % for cross-sectional area
      dsum=10;
      dsum1=10;
      n=-1;
      while ((dsum>1e-6) | (dsum1>1e-6))    % Accuracy set to 0.000001
          sum1 = the_sum;
          sum11= the_sum1;
            n  = n+1;
        d(1,1) = rhor*xs2^2*hspher_1(n, x);
        d(1,2) = ((2*n*(n+1))-xs2^2)*jspher(n,xd2)-4*xd2*der_jspher(n,xd2);
        d(1,3) = 2*n*(n+1)*(xs2*der_jspher(n,xs2)-jspher(n, xs2));
        d(2,1) = -x*der_hspher_1(n,x);
        d(2,2) = xd2*der_jspher(n,xd2);
        d(2,3) = n*(n+1)*jspher(n,xs2);
        d(3,1) = 0;
        d(3,2) = 2*(jspher(n,xd2)-xd2*der_jspher(n,xd2));
        d(3,3) = 2*xs2*der_jspher(n,xs2)+((xs2)^2-2*n*(n+1)+2)*jspher(n,xs2);
      %     
          A(1) = -rhor*(xs2)^2*jspher(n,x);
          A(2) = x*der_jspher(n, x);
          A(3) = 0;
      %
             b = d;
      b(1:3,1) = [A(1) A(2) 0]';
      %
       the_sum = the_sum+((-1)^n)*(2*n+1)*det(b)/det(d);
          dsum = abs(the_sum-sum1);        % Estimate error 
       the_sum1= the_sum1+(2*n+1)*real(det(b)/det(d));
          dsum1= abs(the_sum1-sum11);      % Estimate error
      end
      %
      f(I)   = abs( (2/(complex(0,1)*x) )* the_sum ) ;
      chi(I) = abs( (-2/x^2)* the_sum1);
    end
  ff(:,index) = f(:);
  xx(:,index) = chi(:);
  end
end
%
%----------------------------------------------------------------------
%---                       Functions it calls                       ---           
%----------------------------------------------------------------------
%
function h1 = der_jspher(n,x)
%
%    h1 = der_jspher(n,x)
%
%    Function that returns the derivative of the spherical bessel 
%    function of first kind of order n, with argument x
%
h1 = (n* jspher(n-1, x) - (n+1)* jspher(n+1,x))/(2*n+1) ;
return
end
%
%-----------------------------------------------------------------------
function h1 = der_hspher_1(n,x)
%
%    h1 = der_hspher_1(n,x)
%
%    Function that returns the derivative of the spherical hankel 
%    function of first kind of order n, with argument x
%
h1 = (n* hspher_1(n-1, x) - (n+1)* hspher_1(n+1,x))/(2*n+1) ;
return
end
%
%------------------------------------------------------------------------
function h1 = hspher_1(n,x)
%
%    h1 = hspher_1(n,x)
%
%    Function that returns the Hankel function of first kind
%    of order n, of the argument x
%
h1 = sqrt(pi/(2*x))*besselh( n+0.5, 1, x);
return
end
%
%------------------------------------------------------------------------
function h1 = jspher(n,x)
%
%    h1 = jspher(n,x)
%
%    Function that returns the spherical bessel function of first kind
%    of order n, of the argument x
%
h1 = sqrt(pi/(2*x))*besselj( n+.5, x);
return 
end
%-----------------------------------------------------------------------
