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
##
## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-07-02

## -*- texinfo -*-
## @deftypefn {Function File} {} corr_data (@var{dim}, @var{disc})
## Create a random signal, and impose
## correlation.
##
## @var{dim} controls the size of the field, and
## @var{disc} controls the size of the disc used to filter the field.
##
## This function requires the image processing package.
##
## Reference: THE VARIOGRAM AND ITS ESTIMATION, H. Omre, 
## Geostatistics for Natural Resources Characterization, Part 1, 
## pp. 107-125, 1984 
##
## @end deftypefn
function c = corr_data(dim, disk)

if nargin < 2
    error('Expected 2 arguments')
end

if any(size(disk) ~= 1)
    error('Disc size should be a single integer')
end

a = randn(dim); 
f = fspecial('disk', disk);
c = filter2(f, a);

c = c + abs(min(c(:)));
c = c./max(c(:));

endfunction

%!demo
% cor = corr_data(100, 30);
% surf(cor) 