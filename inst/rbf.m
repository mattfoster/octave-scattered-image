## Copyright (C) 2008 Matt Foster
##  
## Permission is hereby granted, free of charge, to any person obtaining a copy of 
## this software and associated documentation files (the "Software"), to deal in 
## the Software without restriction, including without limitation the rights to 
## use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
## the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##  
## The above copyright notice and this permission notice shall be included in all 
## copies or substantial portions of the Software.
##  
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
## FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
## COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
## IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
## CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


## -*- texinfo -*-
## @deftypefn {Function File} {@var{zz}=} rbf (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
## @deftypefnx {Function File} {@var{zz}=} rbf (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y}, @var{func} ) 
## @deftypefnx {Function File} {@var{zz}=} rbf (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y}, @var{func}, @var{lambda} ) 
## RBF Perform radial basis function interpolation 
##
## Interpolate the scattered values @var{xi}, @var{yi}, @var{zi} 
## at @var{x}, @var{y} using the RBF method.
##
## Functionality is similar to griddata.
##
## Available basis functions (@var{func}):
## @itemize @bullet 
## @item 
## @code{euclidean}:											
## Default. Also fits first order polynomial surface.
##
## @item 
## @code{thin_plate_spline}/@code{biharmonic}:			 
## These are the same. Also fits polynomial surface.
## 
## @item
## @code{gaussian}:
## Use a Gaussian RBF. Constant not implemented.
## 
## @item
## @code{multiquadratic}:
## Use a multiquadratic basis function. Constant not implemented. 
##
## @item
## @code{triharmonic}:
## 
## @end itemize
##
## Regularisation @var{lambda}:
## @itemize @bullet
## @item
## Regularisation parameter see [2]. With regularisation, the interpolation
## condition is relaxed, and the fitted surface will be smoothed. Setting this
## to zero implies interpolation, and high values reduce to a least-squares
## affine model.
## @end itemize 
##
## Output Arguments:
## If multiple output arguments are detected, the number will determine what is returned.
## @enumerate
## @item @var{zz}: Return interpolated output.
## @item @var{zz}, @var{time}: Return time, and interpolated output.
## @item @var{zz}, @var{yy}, @var{zz}: Return all coordinates.
## @item @var{zz}, @var{yy}, @var{zz}, @var{time}: Return all coordinates and time.
## @end enumerate
##
## Example:
## @group
## @example
## [yy, xx] = meshgrid(1:100, 1:100);
## in = corr_data(100, 15);
## aa = rand(100) > 0.95;
## sampled = in.*aa;
## [xi, yi, zi] = find(sampled);
##
## [zz, time] = rbf(xi, yi, zi, xx, yy);
## [zz_tps, time] = rbf(xi, yi, zi, xx, yy, 'thin_plate_spline')
##
## rmse_zz = sqrt(mean((in(:) - zz(:)).^2)) 
## rmse_zz_tps = sqrt(mean((in(:) - zz_tps(:)).^2))
## @end example
## @end group
##
## References: 
## @itemize @bullet	
## @item [1] Surface interpolation with radial basis functions, Carr et al, IEEE
## Trans. on Medical Imaging. Vol 16. No. 1, Feb. 1997
##
## @item [2] Shape matching and object recognition using contexts. Belongie et al,
## IEEE Trans. Pattern Analysis and Machine Intelligence, Vol. 24, No. 24,
## April 2002
## @end itemize
## @seealso{griddata, kriging}
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-26

function [varargout] = rbf(xi, yi, zi, xx, yy, basis_func, lambda, log_fudge)

# Start timer
tic;

# Check input and output argument numbers
error(nargchk(5, 8, nargin));
error(nargchk(1, 4, nargout))

# Change this to control the behaviour of log(0). 
if nargin < 8
  log_fudge = 0;
endif

if nargin >= 6
  basis_func = str2func(basis_func);
else
  basis_func = @euclidean;
endif

[x1,x2]      = meshgrid(yi);
[y1,y2]      = meshgrid(xi);

[basis, poly_order] = basis_func(x1, x2, y1, y2, log_fudge);

if nargin >= 7
  scale = euclidean(x1, x2, y1, y2);
  scale = mean(scale(:));
  lambda = scale.^2 .* lambda;
else
  lambda = 0;
endif

# Perform regularisation (or not)
basis = basis + lambda * eye(size(basis));

# So far I've only seen first order polynomials, or none
if poly_order == 0
  mat = basis;
  ff  = zi;
elseif poly_order == 1 
  q   = [ones(length(xi), 1), xi, yi];
  mat = [basis, q; q', zeros(3)];
  ff  = [zi; zeros(3,1)];
endif

# Things below could get quite memory intensive.
clear basis;

# Get the coefficients.
lam_c = mat \ ff;

# Extract the poly coefficients.
lam   = lam_c(1:end-poly_order*2-1*poly_order);
c     = lam_c(end-poly_order*2:end);
clear lam_c;
clear mat;

# Now reconstruct everything:
if poly_order > 0
  poly = c(1) +  xx.*c(2) + c(3).*yy;
else
  poly = zeros(size(xx));
endif

# Create matric for output
zz = zeros(size(xx));

# This is a vectorised version of the above:
for ii = 1:length(xx(:))

  nm = basis_func( ...
    repmat(xx(ii), [size(xi(:), 1), 1]), ...
    xi(:), ... 
    repmat(yy(ii), [size(yi(:), 1), 1]), ...
    yi(:), ...
    log_fudge);

  zz(ii) = sum(lam .* nm);
end

# Add the polynomial
zz = zz + poly;

# stop timer
time = toc;

# 4 output args -> x, y, z, time
# 3 output args -> x, y, z
# 2 output args -> z, time
# 1 output arg  -> z
if nargout == 4
  varargout{1} = xx;
  varargout{2} = yy;
  varargout{3} = zz;
  varargout{4} = time;
elseif nargout == 3
  varargout{1} = xx;
  varargout{2} = yy;
  varargout{3} = zz;
elseif nargout == 2
  varargout{1} = zz;
  varargout{2} = time;
elseif nargout == 1 
  varargout{1} = zz;
endif

endfunction # rbf 

# Basis functions
# mlint says these are never used.. it's wrong!


# Poly is an integer detailing the order of polynomial to use.

# Euclidean distance -- Linear basis function
function [bf, poly] = euclidean(x1, x2, y1, y2, varargin) ##ok
  bf = sqrt((x1 - x2).^2 + (y1 - y2).^2);
  poly = 1;
endfunction

# Thin plate spline 
function [bf, poly] = thin_plate_spline(x1, x2, y1, y2, varargin) ##ok
  rr = euclidean(x1, x2, y1, y2);
  # Not 100% sure what base to use for the log.
  bf = rr.^2 .* log(rr);
  bf(rr == 0) = varargin{1};
  poly = 1;
endfunction

# Gaussian
function [bf, poly] = gaussian(x1, x2, y1, y2, varargin) ##ok
  rr = euclidean(x1, x2, y1, y2);
  if nargin < 5
    aa = 1;
  else
    aa = varargin{1};
  endif
  bf = exp(aa.*rr.^2);
  # No polynomial here.
  poly = 0;
endfunction

# Multiquadratic
function [bf, poly] = multiquadratic(x1, x2, y1, y2, varargin) ##ok
  rr = euclidean(x1, x2, y1, y2);
  if nargin < 5
    cc = 1;
  else
    cc = varargin{1};
  end
  bf = sqrt(rr.^2 + cc.^2);
  poly = 0;
endfunction

# Triharmonic Spline -- C2 Continuity
function [bf, poly] = triharmonic_spline(x1, x2, y1, y2, varargin) ##ok
  rr = euclidean(x1, x2, y1, y2);
  bf = rr.^4 .* log(rr);
  bf(rr == 0) = varargin{1};
  poly = 0;
endfunction

