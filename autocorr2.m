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
## @deftypefn {Function File} {@var{rho}=} = autocorr2 (@var{x})
## Calculate the 2D Spatial autocorrelation of @var{x} for all lags.
## Output will be a matrix of dimension size(@var{x})).
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-26

function rho = autocorr2(x)
	
	[M,N] = size(x);
	rho   = zeros(M,N);

	x = x - nanmean(x(:)); # Subtract the mean here as part of the normalisation
	v = nanstd(x(:), 1).^2;

	for dx=0:M-1
		for dy=0:N-1
			# Calculate premultiplication 
			pre = 1./((M-dx) .* (N-dy) .* v);	
			# Calculate lag vectors
			lag_x=(0:M-1-dx)+1; 
			lag_y=(0:N-1-dy)+1;
			rho(dx+1,dy+1) = pre .* sum(sum(x(lag_x,lag_y) .* x(lag_x+dx,lag_y+dy)));
		endfor
	endfor
	
endfunction
