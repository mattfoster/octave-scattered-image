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
## @deftypefn {Function File} {@var{ret}=} strat (@var{data}, @var{rate})
## Sample 2-D input @var{data} using stratified sampling, with blocks of size @var{sz}. 
## @var{ret} is an nx3 matrix, corresponing to @code{[row, col, val]}.
##
## Reference:
## @itemize
## @item Spatial Statistics, B. D. Ripley, Wiley Interscience, 2004.
## @end itemize
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-07-04

function [ ret ] = strat ( data, sz )

	if nargin < 2
		sz = 8;
	endif

	samp  = floor(size(data) ./ sz)-1;
	[a,b] = meshgrid([1:samp(1)]*sz, [1:samp(2)]*sz');
	indx  = a + round(rand(samp(1),samp(2)) * sz);
	indy  = b + round(rand(samp(2),samp(1)) * sz);

	indx = indx(:);
	indy = indy(:);

	ret = [indx, indy, fa(sub2ind(size(fa), indx, indy))];

endfunction

## demo
%! a = -0.3 + 4*randn(120);
%! f = fspecial('disk', 15);
%! fa = filter2(f, a);
%! samples = strat(fa);