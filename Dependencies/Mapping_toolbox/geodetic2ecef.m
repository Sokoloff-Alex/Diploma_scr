function [x, y, z] = geodetic2ecef(phi, lambda, h, ellipsoid)
%GEODETIC2ECEF Convert geodetic to geocentric (ECEF) coordinates
%
%   The GEODETIC2ECEF function will be removed in a future release. Use
%   the spheroid/geodetic2ecef method instead.
%
%   [X, Y, Z] = GEODETIC2ECEF(PHI, LAMBDA, H, ELLIPSOID) converts geodetic
%   point locations specified by the coordinate arrays PHI (geodetic
%   latitude in radians), LAMBDA (longitude in radians), and H (ellipsoidal
%   height) to geocentric Cartesian coordinates X, Y, and Z.  The geodetic
%   coordinates refer to the reference ellipsoid specified by ELLIPSOID, an
%   object with SemimajorAxis and Eccentricity properties, or a row vector
%   of the form [semimajor_axis, eccentricity].  H must use the same units
%   as the semimajor axis; X, Y, and Z will be expressed in these units
%   also.
%
%   The geocentric Cartesian coordinate system is fixed with respect to the
%   Earth, with its origin at the center of the ellipsoid and its X-, Y-,
%   and Z-axes intersecting the surface at the following points:
%
%                PHI  LAMBDA
%      X-axis:    0     0      (Equator at the Prime Meridian)
%      Y-axis:    0   pi/2     (Equator at 90-degrees East)
%      Z-axis:  pi/2    0      (North Pole)
%
%   A common synonym is Earth-Centered, Earth-Fixed coordinates, or ECEF.

% Copyright 2004-2012 The MathWorks, Inc.

if isnumeric(ellipsoid)
    % Assume that the input is an ellipsoid vector and extract semimajor
    % axis and flattening.
    a  = ellipsoid(1);
    e2 = ellipsoid(2)^2;
    f = e2 / (1 + sqrt(1 - e2));
else
    a = ellipsoid.SemimajorAxis;
    f = ellipsoid.Flattening;
end

% This function requires input in radians.
inDegrees = false;

[rho, z] = map.geodesy.internal.geodetic2cylindrical(phi, h, a, f, inDegrees);
[x, y] = pol2cart(lambda,rho);
