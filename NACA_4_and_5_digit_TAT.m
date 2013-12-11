%%
% NACA_4_and_5_digit_TAT(alpha, designation) solves thin airfoil theory for NACA
% 4/5 digit airfoils, except for 5 digit reflexed camber lines,
% where 'alpha' is angle of attack in degrees and 'designation' is a 4/5-digit 
% airfoil string.
%
% The five digit camber line is based on curve fits of the constants, 
% so there may be minor discrepancies.
%
% You must enter the airfoil designation as a string,
% if you enter an integer value like 0012, the 00 gets removed.
%
% Example
%   NACA_4_and_5_digit_TAT(5, '2412')
%   will solve for an NACA 2412 at 5 degrees AoA.
%
% Example
%   NACA_4_and_5_digit_TAT(5, '2012')
%   Due to the camber line definitions, this will print a warning and
%   assume that you wanted a NACA 0012 at 5 degrees AoA.
%
% Everything returned wrt alpha is in degrees
%
% Created by: Joseph Derlaga, 11/08/12

function [alpha0L, A0, A1, A2, cm_cx4, cm_LE, x_cp, dcldalpha] = ...
    NACA_4_and_5_digit_TAT(alpha, NACA)

% Uses x (implies x/c) and theta as symbolic values
syms x theta

if numel(NACA) == 4
%% NACA 4 digit mptt
% m  = % max camber
% p  = location of max camber in tenths of chord
% tt = % thickness
% Since this is thin airfoil theory, the thickness doesn't matter
  
  m = str2double(NACA(1))/100;
  p = str2double(NACA(2))/10;

% force symmetric case if p == 0
  if p == 0
    dzdx_fore = 0*(x);
    dzdx_aft  = 0*(x);
    if m ~= 0
      fprintf(strcat('Warning: Assuming you meant 00',NACA(3),NACA(4),'\n'));
    end
  else
    dzdx_fore = (2*m/p^2)*( p - (x) );
    dzdx_aft  = (2*m/(1-p)^2)*( p - (x) );
  end
  
  x_break = p;
else
%% NACA 5 digit mpptt
% l  = design cl / 0.15
% p  = location of max camber in tenths of chord * 2, 
% r  = 1 if reflexed, 0 if not.
% tt = % thickness
% Since this is thin airfoil theory, the thickness doesn't matter

  l = str2double(NACA(1));
  p = str2double(NACA(2))/10;
  r = str2double(NACA(3));
  
  cli = 0.15*l;
  xf = p/2;
  
  m = fzero( @(m)m*(1-sqrt(m/3))-xf, 0.2 );
  q = (3*m - 7*m^2 + 8*m^3 - 4*m^4)/sqrt(m-m^2) ...
    - 1.5*(1-2*m)*(pi/2 - asin(1-2*m));
  
  k1 = 6*cli/q;
  k2 = ( k1*3*(m-xf)^2 - m^3 ) / (1-m)^3;
  
% force symmetric case if p == 0
  if p == 0
    dzdx_fore = 0*(x);
    dzdx_aft  = 0*(x);
    if m ~= 0
      fprintf(strcat('Warning: Assuming you meant 000',NACA(4),NACA(5),'\n'));
      fprintf(strcat('Which is really just a 00',NACA(4),NACA(5),'\n'));
    end
  elseif r == 0
    dzdx_fore =  (k1/6)*( 3*(x)^2 - 6*m*(x) + (3-m)*m^2 );
    dzdx_aft  = -(k1/6)*m^3;
  else
    dzdx_fore = (k1/6)*( 3*((x)-m)^2 - k2/k1*(1-m)^3 - m^3 );
    dzdx_aft  = (k1/6)*( k2/k1*3*((x)-m)^2 - k2/k1*(1-m)^3 - m^3 );
  end

  x_break = m;
end

%% Perform TAT integrations
% Integration substitution, x = (1/2)*(1-cos(theta));
dzdx_fore = subs(dzdx_fore, x, (1/2)*(1-cos(theta)));
dzdx_aft  = subs(dzdx_aft,  x, (1/2)*(1-cos(theta)));

theta_break = acos(1-2*x_break);

% Evaluate the integrals to get the coefficients
alpha0Lrad = eval( (-1/pi) ...
           * ( int(dzdx_fore*(cos(theta)-1), 0, theta_break) ...
             + int(dzdx_aft *(cos(theta)-1), theta_break, pi) ) );

A1 = eval( (2/pi)*( int(dzdx_fore*cos(theta), 0, theta_break) ...
                  + int(dzdx_aft *cos(theta), theta_break, pi) ) );

A2 = eval( (2/pi)*( int(dzdx_fore*cos(2*theta), 0, theta_break) ...
                  + int(dzdx_aft *cos(2*theta), theta_break, pi) ) );

%% Derived coefficients
%alpha0L = radtodeg(alpha0Lrad);
%cl = 2*pi*(degtorad(alpha) - alpha0Lrad);
alpha0L = alpha0Lrad*180.0/pi;
cl = 2*pi*(alpha*pi/180.0 - alpha0Lrad);
A0 = 0.5*(cl/pi-A1);
cm_cx4 = (pi/4)*(A2-A1);
cm_LE  = -(pi/2)*(A0+A1-A2/2);
x_cp   = 0.25*(1+(A1-A2)*pi/cl);

dcldalpha = cl/(alpha-alpha0L);
