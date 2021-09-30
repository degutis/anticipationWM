function [ y ] = lfsr( feedback, start, N_points, f_decimate )
%[ y ] = lfsr( feedback, start, N_points, f_decimate )
%Generate binary sequence from a 32-bit linear feedback shift register (LFSR).
%Inputs:
% feedback: feedback term (from feedback polynomial) which will be XORed
%   with LFSR. Examples for maximal length sequences can be found at
%   http://www.ece.cmu.edu/~koopman/lfsr/index.html (but note these values
%   must be converted from hex to decimal).
% start: starting value for lfsr state. 1 is usually good
% N_points: number of output points to generate
% f_decimate: resample sequence by this factor. 1 is unchanged, 2 takes
%   every second member, etc. Useful for Kasami sequences.
%Note that the mex file is approx 100 times faster.

%    Copyright Travis Wiens 2009
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    If you would like to request that this software be licensed under a less
%    restrictive license (i.e. for commercial closed-source use) please
%    contact Travis at travis.mlfx@nutaksas.com

if nargin<1
    feedback=142;%this would generate an 8-bit maximal length sequence
end
if nargin<2
    start=1;
end
if nargin<3
    N_points=1;
end
if nargin<4
    f_decimate=1;
end

lfsr_state=uint32(start);%get starting value. 
%Note that you may change this and every other uint32 to uint8, uint16 or 
%uint64, depending on your needs
y=zeros(1,N_points);%allocate output

for i=1:N_points
    y(i)=bitand(lfsr_state, uint32(1));%save lsb
    for j=1:f_decimate
        if bitand(lfsr_state, uint32(1))
            lfsr_state=bitxor(bitshift(lfsr_state,-1),uint32(feedback));
        else
            lfsr_state=bitshift(lfsr_state,-1);
        end
    end
end