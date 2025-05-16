function [a,b,c,d]  = plane_3points(A,B,C)

% plane_3points.m calculates the coefficients of the general equation of a
% plane passing through three non-aligned points A, B, and C. The equation
% of the plane is;
% 
% a*x + b*y + c*z + d = 0                                           (eq. 1)
% 
% The coefficients are calculated using the following condition:
% 
%                                                                   (eq. 2)
%     | x-xA,  y-yA,  z-zA|
% det |xB-xA, yB-yA, zB-zA| = 0
%     |xC-xA, yC-yA, zC-zA|
% 
% The coefficients a, b, and c represents the vector normal to the plane.
% For this reason the output is given normalizing the coefficients a, b, c,
% and d by the norm of the vector [a,b,c]'.
% -------------------------------------------------------------------------
% INPUT
% A: x, y, and z coordinates of the first point
% B: x, y, and z coordinates of the second point
% C: x, y, and z coordinates of the third point
% -------------------------------------------------------------------------
% OUTPUT
% a: First coefficient of (eq. 1)
% b: Second coefficient of (eq. 1)
% c: Third coefficient of (eq. 1)
% d: Fourth coefficient of (eq. 1)
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 29 Oct 2016
% California Institute of Technology
% Geological and Planetary Sciences Division

% Solve (eq. 1)
a = (B(2)-A(2))*(C(3)-A(3)) - (B(3)-A(3))*(C(2)-A(2));
b = (B(3)-A(3))*(C(1)-A(1)) - (B(1)-A(1))*(C(3)-A(3));
c = (B(1)-A(1))*(C(2)-A(2)) - (B(2)-A(2))*(C(1)-A(1));
d = -a*A(1)-b*A(2)-c*A(3);
% Norm of the normal vector
n = norm([a,b,c]);
% Normalize the coefficients
a = a/n;
b = b/n;
c = c/n;
d = d/n;