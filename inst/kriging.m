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
## @deftypefn {Function File} {@var{z}=} kriging (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
## @deftypefnx {Function File} {[@var{z}, @var{z_var}]=} kriging (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
## @deftypefnx {Function File} {[@var{x}, @var{y}, @var{z}]=} kriging (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
## @deftypefnx {Function File} {[@var{x}, @var{y}, @var{z}, @var{z_var}]=} kriging (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
## @deftypefnx {Function File} {[@var{x}, @var{y}, @var{z}, @var{z_var}, @var{time}]=} kriging (@var{xi}, @var{yi}, @var{zi}, @var{x}, @var{y} ) 
##
## Interpolate the scattered values @var{xi}, @var{yi}, @var{zi} 
## at @var{x}, @var{y} using the kriging method.
##
## Functionality is similar to griddata.
##
## Code loosely based on recipe from MATLAB recipes for Earth sciences. Martin
## H. Trauth, Springer 2006
##
## @end deftypefn

function [zz, varargout] = kriging(xi, yi, zi, x, y)

	error(nargchk(5, 5, nargin));

	warning off all;

	tic;

	% Create Variogram -- warning.. even this is slow!
	[x1,x2] = meshgrid(xi);
	[y1,y2] = meshgrid(yi);
	[z1,z2] = meshgrid(zi);

	% Calculate all lags.
	D = sqrt((x1 - x2).^2 + (y1 - y2).^2);
	D2 = D.*(diag(xi*NaN)+1);

  % Calculate Semivaiances.
	G = 0.5*(z1 - z2).^2; 

	lag = mean(min(D2));
	hmd = max(D(:))/2;
	max_lags = floor(hmd/lag);
	LAGS = ceil(D/lag);

	for i = 1:max_lags
		SEL = (LAGS == i);
		DE(i) = mean(mean(D(SEL)));
		PN(i) = sum(sum(SEL == 1))/2;
		GE(i) = mean(mean(G(SEL)));
	endfor
	lags=0:max(DE);

	% Fit a spherical variogram model
	mod_func = fit_spherical(DE, GE, var(zi(:)));
	model = mod_func(D);

	% add zeros as the last row and col of the variogram model.
	n = length(xi);
	model(:,n+1) = 1;
	model(n+1,:) = 1;
	% add a zero to the bottom right corner
	model(n+1,n+1) = 0;

	% Invert the variogram model
	model_inv = inv(model);

  % output matrices
	zz = nan(size(x));
	z_var = nan(size(x));

	% Perform the Kriging
	for k = 1:length(zz(:));
		DOR = ((xi - x(k)).^2+(yi - y(k)).^2).^0.5;
		G_R = mod_func(DOR);
		G_R(n+1) = 1; 
		E = model_inv*G_R; 
		zz(k) = sum(E(1:n,1).*zi); 
		z_var(k) = sum(E(1:n,1).*G_R(1:n,1))+E(n+1,1); 
	endfor

	time = toc;

	if nargout == 2 || nargout == 3
		varargout{1} = z_var;
	endif

	if nargout == 3
		varargout{2} = time;
	endif

	if nargout >= 4
		tmp = zz;
		zz = xx;
		varargout{1} = yy;
		varargout{2} = tmp;
		varargout{3} = z_var;
	endif

	if nargout == 5
		varargout{4} = time;
	endif

	warning on all;

endfunction

function mod_func = fit_spherical(DE, GE, var_z)

	% Only fit to the section below 90% of the variance
	ind = 1:min(find(GE > var_z * 0.9));
	coef = [DE(ind);DE(ind).^3]' \ GE(ind)';
	c = sqrt((-0.5 .* (1.5 ./ coef(1)).^-3)./coef(2));
	a = 1.5 .* c / coef(1);

	if isreal(a)
		% Return an anonymous function that models the variogram
		mod_func = @(h) (coef(1) .* h + coef(2) .* h.^3) .* (h <= a) ...
		+ c.*ones(size(h)) .* (h > a);
	else
		mod_func = @(h) (coef(1) .* h + coef(2) .* h.^3);
	endif

endfunction

