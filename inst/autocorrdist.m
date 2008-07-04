## Copyright (C) 2008 Matt Foster
## 
## This program is free software; redistribution and use in source and
## binary forms, with or without modification, are permitted provided that
## the following conditions are met:
## 
##    1.Redistributions of source code must retain the above copyright
##      notice, this list of conditions and the following disclaimer.
##    2.Redistributions in binary form must reproduce the above copyright
##      notice, this list of conditions and the following disclaimer in the
##      documentation and/or other materials provided with the distribution.
## 
## THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
## OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
## HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
## OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
## SUCH DAMAGE.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{lag}=} autocorrdist (@var{xx}, @var{thresh}, @var{smoothing})
## Calculate the lag at whch the 2-D autocorrelation of @var{xx} drops below 
## the threshold @var{thresh}.
## @var{smoothing} sets the size of moving average filter to use for smoothing. 
## Default is 10.
##
## @example
## @group
## [xx,yy] = meshgrid(-pi:0.5:pi);
## rr = sqrt(xx.^2 + yy.^2);
## xx = (2.*cos(r)+rand(size(r)));
## lag = autocorrdist(xx, 1/exp(1));
## @end group
## @end example
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-26

function [ lag ] = autocorrdist (xx, thresh, smoothing)

  error(nargchk(2, 3, nargin));

  if nargin < 3
		smoothing = 10;
	endif

	# First, get the acf. 
	acf = autocorr2(xx);

  # Calculate the lags.
	[rr, cc] = meshgrid(1:size(acf,1), 1:size(acf,2));
	lags = sqrt(rr.^2 + cc.^2);

	[sorted_lags, sort_order] = sort(lags(:));

	smoothed_sorted_acf = filtfilt(ones(smoothing,1)./smoothing, 1, acf(sort_order));
	
	# Find the lowest point that the smoothed acf drops below threshold.
	pos = min(find(smoothed_sorted_acf < thresh));
	lag = sorted_lags(pos);

	# Plot things:
	# plot(sorted_lags, acf(sort_order), 'rx');	
	# hold on,
	# plot(sorted_lags, smoothed_sorted_acf);
	
endfunction
