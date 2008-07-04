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
## @deftypefn {Function File} {@var{v}, @var{range} =} variogram2 (@var{xx}, @var{yy}, @var{zz}, @var{angle_range}, @var{disp_range} )
## Create a variogram from 2D data.
##       @var{xx}, @var{yy}, @var{zz} are coordinates and values. 
##       @var{ang_range} is the range of angles to consider.
##       @var{disp_range} is the size of lag to consider one value.
##       @var{m_ang} is a boolean variable which sets whether angles are modulo pi.
##       @var{alg} is the algorithm to use to construct the variogram, and can be:
## @group
##  'c'  for classis estimation (default)
##  'r' for robust mean based estimation.
##  'm' for robust median based estimation.
##  'n' to use the next argument as a function pointer.
## @end group
##
## See Statistics for Spatial Data, Noel A. C. Cressie, 1991. Wiley and Sons,
## Inc, pages, 69, 75, 215--216.
##
## @end deftypefn

function [v, dist] = variogram2(xx, yy, zz, ang_range, disp_range, m_ang, alg, func)


	if size(ang_range) ~= [1, 2]
		error('Angle range should be 1 x 2');
	end

	if nargin < 4    ang_range = [-pi/4, pi/4];
		warning('Setting default range for angles');
	end

	if nargin < 5
		disp_range = 8;
		warning('Setting default range for displacement');
	end

	if nargin < 7
		estimator = @classical;
	else
		switch alg
		case 'c'
			estimator = @classical;
		case 'r'
			estimator = @robust_mean;
		case 'm'
			estimator = @robust_median;
		case 'n'
			if isa(func, 'function_handle') & nargin == 8
				estimator = func;
			else
				warning('Input argument not a function handle. Using default.');
			end
		otherwise
			warning('Misunderstood estimator algorithm. Using default.');
			estimator = @classical;
		end
	end

	# pairs = nchoosek(1:length(xx), 2);
	# p1 = [xx(pairs(:,1)), yy(pairs(:,1))];
	# p2 = [xx(pairs(:,2)), yy(pairs(:,2))];

	[xp, yp] = meshgrid(1:length(xx));
	pairs = [xp(:), yp(:)];
	p1 = [xx(pairs(:,1)), yy(pairs(:,1))];
	p2 = [xx(pairs(:,2)), yy(pairs(:,2))];

	# Difference each pair to extract lag vectors:
	diff = p1 - p2;
	xd = diff(:,1);
	yd = diff(:,2);
	# Find distances:
	[th, r] = cart2pol(xd, yd);

	# Threshold angles:
	m_ang = 1;
	if m_ang
		th = mod(th, pi);
	end
	th_ind = find(th > ang_range(1) & th < ang_range(2));

	# Threshold displacements
	th_r = [];
	step = disp_range;
	range = 1:ceil(max(r)./step);
	for h = range
		th_r = find(r >= h*(step-1) & r < h*step);
		ind{h} = intersect(th_ind, th_r);
	end

	for ii = 1:length(ind)
		nh = zz(pairs(ind{ii},2)) - zz(pairs(ind{ii},1));
		dist(ii) = mean(r(ind{ii}));
		if all(size(nh)) > 0
			v(ii) = estimator(nh);
		else
			v(ii) = 0;
		end
	end
	
endfunction
# Subfunctions below are various variogram estimators:

# Classical estimator (Cressie, page: 69)
function g = classical(nh)
	g = mean(nh.^2);
endfunction

# Variance based estimator (Bad)
function g = var_est(nh)
	g = var(nh);
endfunction

# Robust mean based estimator (Cressie, 75)
function g = robust_mean(nh)
	g = mean(sqrt(abs(nh))).^4./(0.457+0.494./length(nh));
endfunction

# Robust median based estimator (Cressie, 75)
function g = robust_median(nh)
	g = median(sqrt(abs(nh))).^4./(0.457);
endfunction

# Demo
%! a = -0.3 + 4*randn(120);
%! f = fspecial('disk', 15);
%! fa = filter2(f, a);
%! op = strat(fa, 8);
%! [vc, rc] = variogram2(op(:,1), op(:,2), op(:,3), [0, pi], 2.5, 1, 'c');
%! [vr, rr] = variogram2(op(:,1), op(:,2), op(:,3), [0, pi], 2.5, 1, 'r');
%! [vm, rm] = variogram2(op(:,1), op(:,2), op(:,3), [0, pi], 2.5, 1, 'm');
