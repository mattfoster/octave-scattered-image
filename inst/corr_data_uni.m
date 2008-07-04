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
## @deftypefn {Function File} {} corr_data_uni (@var{dim}, @var{disc})
## Create a random signal, and impose
## correlation. Then remove multinormality by choosing randomly from 10 neighest
## values in a 5x5 neighbourhood.
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

## Author: Matt Foster <matt.p.foster@gmail.com>
## Created: 2008-07-02

function c = corr_data_uni(dim, disk, pad, nh)

% Defaults, these will be overriden if arguments are supplied.

if nargin < 2
    error('Expected at least 2 arguments')
end

if nargin < 3
    pad = 50;
end

if nargin < 4
    nh = 2; % 2*nh+1 square neighbourhood
end

if any(size(disk) ~= 1)
    error('Disc size should be a single integer')
end


% Get some multinormal data
d = corr_data(dim+(2*pad), disk);

% remove multinormality
% Look at local neighbourhood 5x5 mask choose randomly from
% top 10

for ii = pad:size(d,1)-pad-1
    for jj = pad:size(d,nh)-pad-1
        nhood = d(ii-nh:ii+nh, jj-nh:jj+nh);
        sv = sort(nhood(:),'descend');
        pos = ceil(10*rand(1));
        new(ii-pad+1,jj-pad+1) = sv(pos);
    end
end

% Make the histogram look more like omre's
c = abs(reallog(new));
