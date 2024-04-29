function [nr_outputs] = nr_sigmoid(varargin)
% Syntax: There are 2 ways to use this function: "manual" or "automated"
% 
% FOR MANUAL USAGE:
% [nr_outputs] = nr_sigmoid()
% Don't input any parameters.
% Specify your parameter manually below in the parameters section.
%
% FOR AUTOMATED USAGE:
% [nr_outputs] = nr_sigmoid([single_values], rmax, c50, b, n, npoints, x_axis_max)
%
% INPUT PARAMETERS:
% single_values =  varargin{1}: If not empty (not [], but some value.. say
%                              [3]) this specifies that you want to use 
%                              nr_sigmoid to compute the output of value 
%                              fed into the sigmoidal function. 
%                              single_values is a vector and can therefore
%                              contain as many single values as you like. 
%                              When used in this manner, the output values,
%                              will be a vector of numbers corresponding to
%                              the function outputs for single_values. If
%                              you don't want to engage "single_values"
%                              mode, set single_values to [].
% 
% rmax = varargin{2}: defines the asymptotic response maximum as x --> inf
%
% c50 = varargin{3}: defines the point at which the function reaches 50% 
%                    of rmax
%
% b = varargin{4}: allows for a vertical "dc" offset (shifts whole 
%                  function vertically by a value of b..
%
% n = varargin{5}: defines the exponential variable "n" in r(c). 
%                  Affects the curve steepness.
%
% npoints = varargin{6}: specifies the number datapoints sampled in 
%                        the output. *NOTE* : zero is included. 
%                        i.e. if you choose n you will get 0-n  
%                        (n+1 datapoints.. output vector will have
%                        length=n+1). (not necessary in single_values mode)
%
% x_axis_max = varargin{7}: specifies the maximum x value in the
%                           output. Data will range from 0-x_axis_max
%                           along the x-axis. (not necessary in 
%                           single_values mode)
%
% FUNCTION OUTPUTS:
% nr_outputs: The function outputs are saved in a struct called nr_outputs 
%             (or whatever other variable you assign to the output).
%             nr_outputs.xvals: will contain the vector of function
%             x value inputs. nr_outputs.yvals will contain a vector of the
%             y value outputs.
%     
% This function was written to create and plot sigmoidal "Naka Rushton"
% functions. It was, however, intentionally written in a generalized
% manner such that, theoretically, different functions can also be used as 
% well by adjusting the "function definition" parameters below.
%
% Written by Ethan Duwell, PhD in April 2022 
% as a Post-doc in Adam Greenberg's Lab
%
%% Parameters


if nargin == 0  % Params within this if statement need to be specified 
                % manually below if no input parameters are fed in...
                % NOTE: IN LIEU OF ANY INPUT PARAMS, THESE WILL ACT LIKE DEFAULTS!!!
single_values =[]; % Leave empty if you do not wish to engage single_values
                   % mode. If not empty, the program will calculate the
                   % sigmoidal function outputs for only the values in this
                   % vector.

plot_figs = 1; % if 1, nr_sigmoid will plot figures 

% Naka Rushton Sigmoid Curve Params.. 
% (for r(c) below under "function definition parameters")
rmax = 0.5;  % defines the asymptotic response maximum as x --> inf
c50 = 1.0; % defines the point at which the function reaches 50% of rmax
b = 0.5;     % allows for a vertical "dc" offset (shifts whole function vertically by a value of b..
n = 2;     % defines the exponential variable "n" in r(c). Affects the curve steepness.

%%%%%%%%%%%%%%%%% Continuous Function Plotting Params %%%%%%%%%%%%%%%%%
cont_min = 0; % minimum value for continuous plots made with "fplot"
cont_max = 2; % maximum value for continuous plots made with "fplot"

%%%%%%%%%%%%%%%%% Quantized Function Input Params %%%%%%%%%%%%%%%%%
npoints = 100; % note: zero is included. IE if you choose, 100 you will get 0-100 (101 points).
x_axis_max = 2; % Specifies the maximum X value.
else
    single_values =  varargin{1};
    rmax = varargin{2};
    c50 = varargin{3};
    b = varargin{4};
    n = varargin{5};
    if nargin > 5
        npoints = varargin{6};
        x_axis_max = varargin{7};
    end
end

% given the number of datapoints sampled (npoints) and the x_axis_max make
% a vector x that is npoints long ranging from 0-x_axis_max
if isempty(single_values)
x=((0:npoints)/npoints)*x_axis_max;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Function Definition Parameters %%%%%%%%%%%%%%%%%
syms c
r(c) = rmax*((c^n)/(c50^n + c^n)) + b; % naka rushton model function
sigmoid = matlabFunction(r(c)); % assign it as a call-able matlab function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code
clear nr_outputs
clear out
% not single_value mode
if isempty(single_values)
    out = zeros(1, npoints);
    for ii = 1:npoints+1
        out(1,ii) = sigmoid(x(1,ii));
    end
    %y = out;
    out = vertcat(x,out);
end

% single_value mode
if ~isempty(single_values)
    out = zeros(1, length(single_values));
    for ii = 1:length(single_values)
        out(1,ii) = sigmoid(single_values(1,ii));
    end
    %y = out;
    out = vertcat(single_values,out);
end

%% Plots

if exist("plot_figs","var") 
if plot_figs == 1

    % plot of quantized datapoints in the out vector
    figure;
    gcf;
    xx = out(1,1:end);
    yy = out(2,1:end);
    plot(xx,yy);
    axis([0 x_axis_max 0 (rmax+b)]);
    
    % make a continuous plot of the function using fplot
    figure;
    gcf;
    fplot(sigmoid,[cont_min cont_max]);

end
end
%% WRAP UP 
% save desired outputs in a data structure called nr_outputs
if exist("nr_outputs","var") == 1
clear nr_outputs
end
nr_outputs = out;

end