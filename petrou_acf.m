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
## @deftypefn {Function File} {@var{c}=} petrou_acf (@var{m})
## Calculate the 2D Spatial Autocorrelation of @var{m} for all lags.
## Output @var{c} will be a matrix prod(size(m))
#
## As described in Image Processing: The Fundamentals 
## by Maria Petrou and Panagiota Bosdogianni, pp 98.
##
## WARNING: Outputs probably don't have the form you expect!
##
## Output will have the form:
## @group
## A B C D E F G H I
## B A B J D E K G H
## C B A L J D M K G
## D J L A B C D E F
## E D J B A B J D E
## F E D C B A L J D
## G K M D J L A B C
## H G K E D J B A B
## I H G F E D C B A
## A = average squared element
## B = average product of vertical neighbours
## C = average product of vertical neighbours once removed
## D = average product of horizontal neighbours
## E = average product of diagonal neighbours
## F = average procuct of values with vertical distance of 2, horizontal 1
## G = average product of values with horizontal removal of 2
## H = average product of values with horizontal removal of 2, vertical 1
## I = average product of diagonal neighbours once removed
## J = average product of off-diagonal neighbours
## K = average product of values with horizontal removal of 2, vertical -1
## L = average product of values with horizontal removal of -1, vertical 2
## M = values with off-diagonal neighbours once removed.
## @end group
## @end deftypefn

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-06-26

function c = petrou_acf(m)

	if nargin < 1
		error('Not enough arguments')
	endif

	N = size(m, 1);

	for i = 1:prod(size(m));
		for j = 1:prod(size(m))
			k0 = mod(i - 1, N) - mod(j - 1, N);
			l0 = (i - j - k0) ./ N;

			s1_l = max(1, 1 + k0);
			s1_u = min(N, N + k0);
			s2_l = max(1, 1 + l0);
			s2_u = min(N, N + l0);

			t = 0;

			for s1i = s1_l:s1_u
				for s2i = s2_l:s2_u
					t = t + m(s1i, s2i) * m(s1i - k0, s2i - l0); 
				endfor
			endfor

			c(i,j) = t./N.^2;
		endfor
	endfor

endfunction

%!m = repmat([1,2,1], 3, 1);
%!test_a = acf(m);
%!correct = [
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    666.666666666667e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    444.444444444444e-003
%!    333.333333333333e-003
%!    222.222222222222e-003
%!    111.111111111111e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    222.222222222222e-003
%!    333.333333333333e-003
%!    222.222222222222e-003
%!    666.666666666667e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!    444.444444444444e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    111.111111111111e-003
%!    222.222222222222e-003
%!    333.333333333333e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    444.444444444444e-003
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    666.666666666667e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    444.444444444444e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    444.444444444444e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    666.666666666667e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!    444.444444444444e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    333.333333333333e-003
%!    222.222222222222e-003
%!    111.111111111111e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    444.444444444444e-003
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    666.666666666667e-003
%!    222.222222222222e-003
%!    333.333333333333e-003
%!    222.222222222222e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!    1.33333333333333e+000
%!    111.111111111111e-003
%!    222.222222222222e-003
%!    333.333333333333e-003
%!    444.444444444444e-003
%!    888.888888888889e-003
%!    1.33333333333333e+000
%!    666.666666666667e-003
%!    1.33333333333333e+000
%!    2.00000000000000e+000
%!];
%!
%!assert(all(correct - test_a(:) < 1e-6));
