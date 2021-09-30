/*
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
 */

/*
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
 */

/* $Revision: 1.10 $ */
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	F_IN	prhs[0]
#define	S_IN	prhs[1]
#define	N_IN	prhs[2]
#define	FD_IN	prhs[3]

/* Output Arguments */

#define	Y_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.14159265



static void lfsr(double *y, unsigned int feedback, unsigned int start, unsigned int N_points, unsigned int f_decimate)
{
    unsigned int lfsr_state;
    unsigned int i,j;
    lfsr_state=start;
    
    for (i=0;i<N_points;i++)
    {
        y[i]=(double) (lfsr_state & 0x01);
        for (j=0;j<f_decimate;j++)
        {
        if (lfsr_state & 1)  { lfsr_state = (lfsr_state >> 1) ^ feedback; }
            else        { lfsr_state = (lfsr_state >> 1);}
        }
        
    }
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
   
    unsigned int feedback;
    unsigned int start;
    unsigned int N_points;
    unsigned int f_decimate;

    double *y;
    
    /* Check for proper number of arguments */
    
    if (!((nrhs == 4)||(nrhs == 3))) { 
	mexErrMsgTxt("Three or Four input arguments required: feedback, start, N_points, f_decimate"); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    


    
    feedback = *mxGetPr(F_IN); 
    start=*mxGetPr(S_IN);
    N_points=*mxGetPr(N_IN);
    if (nrhs == 3) {
        f_decimate=1;
    } else {
        f_decimate=*mxGetPr(FD_IN);
    }      
    /* Create a matrix for the return argument */ 
    Y_OUT = mxCreateDoubleMatrix(1, N_points, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    y = mxGetPr(Y_OUT);
    /* Do the actual computations in a subroutine */
    lfsr(y,feedback,start,N_points,f_decimate); 
    return;
    
}


