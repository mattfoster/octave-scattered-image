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
## @deftypefn {Function File} {@var{c}=} petrou_acf_alt (@var{m})
## Calculate the 2D Spatial Autocorrelation of @var{m} for all lags.
## Output @var{c} will be a matrix prod(size(m))
##
## As described in Image Processing: The Fundamentals 
## by Maria Petrou and Panagiota Bosdogianni
##
## WARNING: Outputs probably don't have the form you expect!
##
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-26

function c = petrou_acf_alt(aa)
	
	if nargin < 1
		error('Expected 1 input argument');
	endif

	[sx,sy] = size(aa);
	orig = aa;

	# Pad input array
	aa = zeros(3*size(orig)-2);
	aa(sx:2*sx-1, sy:2*sy-1) = orig; 

	i_lim = sx-1;
	j_lim = sy-1;

	# Calculate covariances for each lag
	for ii=-i_lim:i_lim
		for jj=-j_lim:j_lim
			l1i = ii+1;
			l1j = jj+1;
			c(ii+i_lim+1,jj+j_lim+1) = ...
			sum(sum(aa(l1i+i_lim:l1i+2*i_lim, l1j+j_lim:l1j+2*j_lim) .*  orig));
		endfor
	endfor

	# Normalise covariances
	cv = c./prod(size(orig));
	clear c;

	# Get lag positions
	ri = cols(orig);
	ci = rows(orig);
	ri = ri(:);
	ci = ci(:);

	# Create final ACF matrix (as a col)
	for ii=1:length(ri)
		c(ii) = cv(ci(ii) , ri(ii));
	endfor

	# Reshape to final matrix
	sz = size(ci); 
	c = reshape(c, sz(1), sz(2));
endfunction

% Functions to create matrices of lags:
function r = cols(aa)
	c = repmat(1:size(aa,1), size(aa,1), 1)
	c = c(:);
	r = lag_mat(c) + size(aa,1);
endfunction

function c = rows(aa)
	a = repmat(1:size(aa,2), 1, size(aa,2))
	c = lag_mat(a) + size(aa,2);
endfunction

function d = lag_mat(in)
	[a,b] = meshgrid(in, in);
	d = a - b;
endfunction
