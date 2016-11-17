function [phi, lambda, h] = ecef2geodetic(x, y, z, ellipsoid)
%ECEF2GEODETIC Convert geocentric (ECEF) to geodetic coordinates
%
%   The ECEF2GEODETIC function will be removed in a future release. Use
%   the spheroid/ecef2geodetic method instead.
%
%   [PHI, LAMBDA, H] = ECEF2GEODETIC(X, Y, Z, ELLIPSOID) converts point
%   locations in geocentric Cartesian coordinates, stored in the
%   coordinate arrays X, Y, Z, to geodetic coordinates PHI (geodetic
%   latitude in radians), LAMBDA (longitude in radians), and H (height
%   above the ellipsoid). The geodetic coordinates refer to the reference
%   ellipsoid specified by ELLIPSOID, an object with SemimajorAxis and
%   Eccentricity properties, or a row vector of the form [semimajor_axis,
%   eccentricity]. X, Y, and Z must use the same units as the semimajor
%   axis;  H will also be expressed in these units.  X, Y, and Z must have
%   the same shape; PHI, LAMBDA, and H will have this shape also.
%
%   For a definition of the geocentric system, also known as
%   Earth-Centered, Earth-Fixed (ECEF), see the help for GEODETIC2ECEF.

% Copyright 2005-2012 The MathWorks, Inc.

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

[lambda, rho] = cart2pol(x,y);
[phi, h] = map.geodesy.internal.cylindrical2geodetic(rho, z, a, f, inDegrees);
