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

## -*- textinfo -*-
## @deftypefn {Function File} {} first_order_nc (@var{si}, @var{cm}, @var{rr})
## Perform first order normalised convolution on an image.
##
## @var{si} is the irregularly sampled image matrix.
## @var{cm} is the confindence map (logical).
## @var{rr} is the radius of the Gaussian filter to use.
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-17


function [ff, dx, dy] = first_order_nc(si, cm, rr)
## function [xx, yy] = first_order_nc(si, cm, rr)
##	Perform first order normalised convolution on an image.
#
##	si is the irregularly sampled image.
##	cm is the confindence map.
## rr is the radius of the Gaussian filter to use.
#
# 	See also NDCMASK, ADAPTIVE_NC

error(nargchk(2, 3, nargin));

if nargin < 3
	rr = 11;
end

[m, m_x, m_y, m_xy, m_xx, m_yy] = generate_filters(rr);

cm 	  = double(cm); 

# Calculate coefficients for the first matrix.
ca_1  = conv2(cm, m, 'same');
ca_x  = conv2(cm, m_x, 'same');
ca_y  = conv2(cm, m_y, 'same');
ca_xx = conv2(cm, m_xx, 'same');
ca_yy = conv2(cm, m_yy, 'same');
ca_xy = conv2(cm, m_xy, 'same');

# Calculate coefficients for the second matrix.
sa_1 = conv2(si, m, 'same');
sa_x = conv2(si, m_x, 'same');
sa_y = conv2(si, m_y, 'same');

ff = zeros(size(si));
dx = zeros(size(si));
dy = zeros(size(si));

for ii = 1:length(si(:))
	# Construct the denominator
	D = [ca_1(ii),	ca_x(ii),		ca_y(ii); ...
			 ca_x(ii),	ca_xx(ii),	ca_xy(ii);...
			 ca_y(ii),	ca_xy(ii),	ca_yy(ii)];

	# Construct the numerator
	N = [sa_1(ii), sa_x(ii), sa_y(ii)]';

  # Calculate the outputs 
	tmp = pinv(D) * N;
	ff(ii) = tmp(1);
	dx(ii) = tmp(2);
	dy(ii) = tmp(3);
endfor

endfunction

# Generate filters for first-order NC
function [m, m_x, m_y, m_xy, m_xx, m_yy] = generate_filters(dim);

	[x, y] = meshgrid(-dim:dim, -dim:dim);	
	
	# dim = 2*sig+1;
	sig = (dim - 1)./2;
	
	m = fspecial('gaussian', [2*dim+1, 2*dim+1], sig);

	m_x = m.*x;
	m_y = m.*y;
	m_xy = m.*(x.*y);
	m_xx = m.*(x.^2);
	m_yy = m.*(y.^2);

endfunction