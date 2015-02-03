function [varargout] = AverageAbsProfiles(averageOver,varargin)
%*************************************************************************
% Function to average a set of ABS profiles. 
% Passed
%       AverageOver
%       Any number of 2D Arrays to average over
% 
% Returns 
%       Same number of arrays as Input
%*************************************************************************

        %Calculte the num of Profiles in the result
    for n=1:size(varargin,2)
        V = varargin{n};
        numProfiles = fix(size(V,2)/averageOver);
        Vout = zeros(size(V,1),numProfiles);
        for N = 1:numProfiles            
            Vout(:,N) = mean(V(:,(N-1)*averageOver +1 :(N-1)*averageOver+averageOver),2);
        end
        varargout{n} = Vout;
    end
    
            
end