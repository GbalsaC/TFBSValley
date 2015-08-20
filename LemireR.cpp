// input: Array a, integer window width w
// output: arrays maxval and minval
// buffer: lists U and L
// requires: STL for deque support

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <deque>
#include <fstream>      // fstream
#include <string>
#include <sstream>
#include <algorithm>    // copy
#include <math.h>
#include <iterator>     // ostream_operator
using namespace Rcpp;
using namespace std;
typedef unsigned int uint;

//Allocate memory for size, window and pointer for dynamic memory allocation
double *maxval;
double *minval;

long sizeArray;
long index, index2;
//std::vector<double> a;
std::deque<int> U, L;
int n;


// [[Rcpp::export]]
NumericMatrix lemireR(NumericVector a, int w, int offset, double Depth){
        index = 0; index2 = 0;
        sizeArray = a.size() + 1;
	//A = new double[sizeArray];
        NumericVector A1(sizeArray);
        NumericVector A2(sizeArray);        
	NumericMatrix A(a.size(), 5);
	minval = new double[sizeArray];
	maxval = new double[sizeArray];



	for (long i = 1; i < a.size(); ++i)
	{
		if (i >= w)
		{
			maxval[i - w] = a[U.size()>0 ? U.front() : i - 1];
			minval[i - w] = a[L.size()>0 ? L.front() : i - 1];
		}
		if (a[i]>a[i - 1])
		{
			L.push_back(i - 1);
			if (i == w + L.front())
			{
				L.pop_front();
			}
			while (U.size()>0)
			{
				if (a[i] <= a[U.back()])
				{
					if (i == w + U.front())
					{
						U.pop_front();
					}
					break;
				}
				U.pop_back();
			}
			//End While
		}
		else
		{
			U.push_back(i - 1);
			if (i == w + U.front())
			{
				U.pop_front();
			}
			while (L.size()>0)
			{
				if (a[i] >= a[L.back()])
				{
					if (i == w + L.front())
					{
						L.pop_front();
					}
					break;
				}//END IF
				L.pop_back();
			}//END While
		}//End If-Else
	}//END FOR
	maxval[a.size() - w] = a[U.size()>0 ? U.front() : a.size() - 1];
	minval[a.size() - w] = a[L.size()>0 ? L.front() : a.size() - 1];

	//myfile is optionally only usable as csv output
	for (long l = 1; l < sizeArray; l++){
		if (l < sizeArray - (w))
		{
        		A(l-1+(w/2), 0) = maxval[l];
                        //A(0, l-1) = maxval[l];
        		//myfile << maxval[l] << ",";
		}
                else if(l < sizeArray - (w/2))
		{
			A(l-1+(w/2), 0) = 0;
			//myfile << 0 << ",";
		}
                /* 
                */
	}
	//myfile << "\n";
        /*
        // Minval filter is muted as output since a "custom" implementation is in order
	for (long k = 1; k < sizeArray; k++){
				//if (k < sizeArray - w)
                if (k < sizeArray - (w))
		{
			A(2, k-1+(w/2)) = minval[k];
            //A(1, k-1) = minval[k];
			//myfile << minval[k] << ",";
		}
                
		else if (k < sizeArray - (w/2))
		{
			A(2, k-1+(w/2)) = 0;
			//myfile << 0 << ",";
		}
                
                 
	}
        */
        n = w + offset;
        // C++ code for V(pos) in Ramsey's Matlab 
	for (long pos = 1; pos < sizeArray-1; pos++){
                if ((pos >= n) && (pos < (sizeArray - w))){
			A(pos, 1) = min(A(lrint(pos + offset + w / 2), 0), A(lrint(pos - offset - w / 2),0));
                        
		}
		else {
			A(pos, 1) = 0;
		}
                //Create logical comparison of A(pos,1) with a[]
                //Consider previous pos if, if 0 then it's false for all comparisons
                // Of the type A(pos,1)*Depth > a[pos]
                if((A(pos, 1)*Depth)>a[pos]){
                        A(pos, 2)=A(pos, 1);
                }
                
                
	}
        

	//myfile.close();
	delete[] maxval;
	delete[] minval;
	return A;
}
	