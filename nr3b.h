#ifndef _NR3_H_
#define _NR3_H_

// NR3B Numerical Recipes class library with augmented functionality
// This file nr3b.h is backwards compatible with nr3.h and nr3a.h

#define _CHECKBOUNDS_ 1
#ifdef _DEBUG //NR3A
//#define _USENRERRORCLASS_ 1 // enable for proper traceback in Win Debug, disable for proper message in Release
#endif

#define _TURNONFPES_ 1


// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <string>
#include <ctype.h>
#include <stdarg.h>
#include <typeinfo>
#include <float.h>
#include <cstdint>
#include <exception>
#include <thread>  // for sleep()
#include <chrono>  // for sleep()


using namespace std;

// macro-like inline functions

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

// exception handling

class THROWexception : public exception {
	string message;
public:
	THROWexception(const string& msg, const char* file, int line) {
		message = msg + " in file " + file + " at line " + to_string(line);
	}
	const char* what() const noexcept override {
		return message.c_str();
	}
};
#define THROW(msg) throw THROWexception(msg, __FILE__, __LINE__)
// redefine the above if you want different error handling
// DNAcode and RSecc call THROW(const char *msg) on all errors

// clock
class NRclock {
public:
	// times in milliseconds
	clock_t first, previous, now;
	NRclock() : first(clock()), previous(first) {}
	double lap() { // call before and after desired code
		now = clock();
		double laptime = 1000.*double(now - previous) / CLOCKS_PER_SEC ;
		previous = now;
		return laptime;
	}
	double elapsed() {return 1000.*double(clock() - first) / CLOCKS_PER_SEC;} 
};

// Vector and Matrix Classes

template <class T> class NRmatrix; // forward declarations
template <class T> class NRvectorview; 
template <class T> class NRmatrixview; 

template <class T>
class NRvector {
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(int n, const T &a);	//initialize to constant value
	NRvector(int n, const T *a);	// initialize to array
	NRvector(int n, char *sep, ...); // initialize to varargs (sep is arbitrary string like "values")
	NRvector(const NRvector &rhs);	// copy constructor
	template <class U> NRvector(const NRvector<U> &rhs);	// conversion operator
	NRvector & operator=(const NRvector<T> &rhs);	//assignment
	NRvector & operator=(const T rhs);	//scalar assignment
	typedef T value_type; // make T available externally
	inline T & operator[](int i);	//i'th element
	inline const T & operator[](int i) const;
	inline int size() const;
	inline int end() const;
	NRvector<T> &resize(int newn, bool preserve=false); // resize 
	NRvector<T> &assign(int newn, const T &a); // resize and assign a constant value
	T maxval() const;
	T maxval(int &loc) const;
	T maxval(int &loc1, int &loc2) const;
	T minval() const;
	T minval(int &loc) const;
	T minval(int &loc1, int &loc2) const;
	T sum() const;
	bool contains(const T val) const;
	int index(const T val) const;
	NRmatrixview<T> asrowmat(); // vector as 1 x n matrix
	const NRmatrix<T> asrowmat() const;
	NRmatrixview<T> ascolmat(); // vector as n x 1 matrix
	const NRmatrix<T> ascolmat() const;
	NRvector<T> cat(const NRvector<T> &b) const;
	NRvector<T> cat(T b) const;
	NRvector<T> gather(const NRvector<int> &locs) const;
	template <class U> NRvector<T> &scatter(const NRvector<int> &locs, const NRvector<U> &vals);
	template <class U> NRvector<T> &scatter(const NRvector<int> &locs, const U val);
	NRvector<int> find() const; // return locs of nonzero entries for use with scatter/gather
	NRvectorview<T> slice(int beg, int end);
	const NRvectorview<T> slice(int beg, int end) const;
	NRvector<T> cumsum() const;
	~NRvector();
	template <class U> friend class NRmatrix;
	template <class U> friend class NRvectorview;
	template <class U> friend class NRmatrixview;
};

// NRvector definitions

template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL) {
	for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL) {
	for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(int n, char *sep, ...) : nn(n), v(n>0 ? new T[n] : NULL) {
    va_list ap;
    int j;
    va_start(ap, sep);
    for (j=0;j<n;j++) v[j] = va_arg(ap, T);
    va_end(ap);
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL) {
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T> template <class U>
NRvector<T>::NRvector(const NRvector<U> &rhs) : nn(rhs.size()), v(nn>0 ? new T[nn] : NULL) {
	for(int i=0; i<nn; i++) v[i] = T(rhs[i]);
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs) {
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
	if (this != &rhs) {
		if (nn != rhs.size()) {
			if (v != NULL) delete [] (v);
			nn=rhs.size();
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const T rhs) {
	for (int i=0;i<nn;i++) v[i] = rhs;
	return *this;
}

template <class T>
inline T & NRvector<T>::operator[](const int i)	{ //subscripting 
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	THROW("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[](const int i) const	{ //subscripting
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	THROW("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRvector<T>::size() const { return nn; }

template <class T>
inline int NRvector<T>::end() const { return nn-1; }

template <class T>
NRvector<T>& NRvector<T>::resize(int newn, bool preserve) {
	if (newn != nn) {
		if (preserve) {
			int i,nmin = MIN(nn,newn);
			T *vsave = v;
			v = newn > 0 ? new T[newn] : NULL;
			for (i=0;i<nmin;i++) v[i] = vsave[i];
			for (i=nmin;i<newn;i++) v[i] = T(0);
			if (vsave != NULL) delete[] (vsave);
			nn = newn;
		} else {
			if (v != NULL) delete[] (v);
			nn = newn;
			v = nn > 0 ? new T[nn] : NULL;
		}
	}
	return *this;
}

template <class T>
NRvector<T>& NRvector<T>::assign(int newn, const T& a) {
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0;i<nn;i++) v[i] = a;
	return *this;
}

template <class T>
T NRvector<T>::minval(int &loc1, int &loc2) const {
	int i;
	T val = numeric_limits<T>::max();
	for (i=0;i<nn;i++) if (v[i] <= val) {
		loc2 = i;
		if (v[i] < val) {
			val = v[i];
			loc1 = i;
		}
	}
	return val;
}
template <class T>
T NRvector<T>::minval(int &loc) const {
	int loc2;
	return minval(loc,loc2);
}
template <class T>
inline T NRvector<T>::minval() const {
	int loc1,loc2;
	return minval(loc1,loc2);
}

template <class T>
T NRvector<T>::maxval(int &loc1, int &loc2) const {
	int i;
	T val = numeric_limits<T>::lowest();
	for (i=0;i<nn;i++) if (v[i] >= val) {
		loc2 = i;
		if (v[i] > val) {
			val = v[i];
			loc1 = i;
		}
	}
	return val;
}
template <class T>
T NRvector<T>::maxval(int &loc) const {
	int loc2;
	return maxval(loc,loc2);
}
template <class T>
inline T NRvector<T>::maxval() const {
	int loc1,loc2;
	return maxval(loc1,loc2);
}

template <class T>
T NRvector<T>::sum() const {
	int i;
	T ans = T(0);
	for (i=0;i<nn;i++) ans += v[i];
	return ans;
}

template <class T>
bool NRvector<T>::contains(const T val) const {
	int i;
	for (i = 0; i < nn; i++) if (v[i] == val) return true;
	return false;
}

template <class T>
int NRvector<T>::index(const T val) const {
	int i;
	for (i = 0; i < nn; i++) if (v[i] == val) return i;
	return -1;
}

template <class T>
const NRmatrix<T> NRvector<T>::asrowmat() const {
	return NRmatrix<T>(1,nn,v); // copy data from v
}

template <class T>
NRmatrixview<T> NRvector<T>::asrowmat() {
	NRmatrixview<T> ans;
	ans.nn = 1;
	ans.mm = nn;
	ans.v = new T*[1];
	ans.v[0] = v;
	return ans;
}

template <class T>
const NRmatrix<T> NRvector<T>::ascolmat() const {
	return NRmatrix<T>(nn,1,v); // copy data from v
}

template <class T>
NRmatrixview<T> NRvector<T>::ascolmat() {
	NRmatrixview<T> ans;
	ans.nn = nn;
	ans.mm = 1;
	ans.v = (nn > 0 ? new T*[nn] : NULL);
	for (int i=0;i<nn;i++) ans.v[i] = v+i;
	return ans;
}

template <class T>
NRvector<T> NRvector<T>::cat(const NRvector<T> &b) const {
	int i, nb = b.size();
	NRvector<T> ans(nn+nb);
	for (i=0;i<nn;i++) ans[i] = v[i];
	for (i=0;i<nb;i++) ans[i+nn] = b[i];
	return ans;
}

template <class T>
NRvector<T> NRvector<T>::cat(T b) const {
	int i;
	NRvector<T> ans(nn+1);
	for (i=0;i<nn;i++) ans[i] = v[i];
	ans[nn] = b;
	return ans;
}

template <class T>
NRvector<T> NRvector<T>::cumsum() const {
	int i;
	NRvector<T> ans(nn);
	if (nn>0) ans[0] = v[0];
	for (i=1;i<nn;i++) ans[i] = ans[i-1] + v[i];
	return ans;
}


template <class T>
NRvector<T> NRvector<T>::gather(const NRvector<int> &locs) const {
	int i, n = locs.size();
	NRvector<T> ans(n);
	for (i=0;i<n;i++) ans[i] = operator[](locs[i]); // so that range checking is done if enabled
	return ans;
}

template <class T> template <class U>
NRvector<T>& NRvector<T>::scatter(const NRvector<int> &locs, const NRvector<U> &vals){
	int i, n = locs.size();
	if (vals.size() != n) THROW("arguments in set must be same length");
	for (i=0;i<n;i++) operator[](locs[i]) = T(vals[i]); // so that range checking is done if enabled
	return *this;
}

template <class T> template <class U>
NRvector<T>& NRvector<T>::scatter(const NRvector<int> &locs, const U val) {
	int i, n = locs.size();
	for (i=0;i<n;i++) operator[](locs[i]) = T(val); // so that range checking is done if enabled
	return *this;
}

template <class T>
NRvector<int> NRvector<T>::find() const {
	int i,n;
	for (n=0,i=0;i<nn;i++) if (v[i]) ++n;
	NRvector<int> ans(n);
	for (n=0,i=0;i<nn;i++) if (v[i]) ans[n++] = i;
	return ans;
}

template <class T>
NRvector<T>::~NRvector() {
	//printf("destr %d\n", nn);
	if (v != NULL) delete[] (v);
}

template <class T>
class NRvectorview : public NRvector<T> {
	using NRvector<T>::nn;
	using NRvector<T>::v;
public:
	NRvectorview();
	NRvectorview(const NRvectorview<T> &rhs);
	NRvectorview<T> & operator=(const NRvectorview<T> &rhs);
	NRvectorview<T> & operator=(const NRvector<T> &rhs);
	NRvectorview<T> & operator=(const T rhs);
	friend class NRvector<T>;
	friend class NRmatrix<T>;
	~NRvectorview(); 
};

template <class T>
NRvectorview<T>::NRvectorview() {}

template <class T>
NRvectorview<T>::NRvectorview(const NRvectorview<T> &rhs) { // make shallow copy explicit
	NRvector<T>::nn = rhs.nn;
	NRvector<T>::v = rhs.v;
};

template <class T>
NRvectorview<T> & NRvectorview<T>::operator=(const NRvectorview<T> &rhs) {
	if (this != &rhs) {
		if (nn != rhs.nn) THROW("incompatible sizes in slice assignment");
		for (int i=0; i<nn; i++) v[i]=rhs[i];
	}
	return *this;
}

template <class T>
NRvectorview<T> & NRvectorview<T>::operator=(const NRvector<T> &rhs) {
	if (nn != rhs.nn) THROW("incompatible sizes in slice assignment");
	for (int i=0; i<nn; i++) v[i]=rhs[i];
	return *this;
}

template <class T>
NRvectorview<T> & NRvectorview<T>::operator=(const T rhs) {
	for (int i=0; i<nn; i++) v[i]=rhs;
	return *this;
}

template <class T>
NRvectorview<T>::~NRvectorview() { NRvector<T>::v = NULL; } // prevent NRvector destructor from deleting underlying data

template <class T>
NRvectorview<T> NRvector<T>::slice(int begoff, int endoff){
	NRvectorview<T> ans;
	ans.nn = MAX(0,endoff - begoff + 1);
	ans.v = v + begoff;
	return ans;
}

// end of NRvector definitions

template <class T>
class NRmatrix {
	int nn;
	int mm;
	T **v;
public:
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(int n, int m, const T &a);	//Initialize to constant
	NRmatrix(int n, int m, const T &a, char *dum);	//Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(int n, int m, char *sep, ...); // Initialize to varags
	NRmatrix(const NRmatrix &rhs);		// copy constructor
	template <class U> NRmatrix(const NRmatrix<U> &rhs);		// conversion operator
	NRmatrix & operator=(const NRmatrix<T> &rhs);	//assignment
	NRmatrix & operator=(const T rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](int i);	//subscripting: pointer to row i
	inline const T* operator[](int i) const;
	inline int nrows() const;
	inline int ncols() const;
	inline int rend() const;
	inline int cend() const;
	NRmatrix<T> &resize(int newn, int newm, bool preserve=false); // resize
	NRmatrix<T> &assign(int newn, int newm, const T &a); // resize and assign a constant value
	T maxval() const;
	T maxval(int &rloc, int &cloc) const;
	T minval() const;
	T minval(int &rloc, int &cloc) const;
	T sum() const;
	NRvectorview<T> asvec(); // can be an l-value
	NRvectorview<T> row(int irow); // row() can be an l-value (this might be too confusing)
	const NRvector<T> row(int irow) const;
	NRvector<T> col(int jcol) const;
	NRvector<T> diag(int jdiag=0) const;
	NRmatrix<T> &setrow(int irow, const NRvector<T> &row);
	NRmatrix<T> &setcol(int jcol, const NRvector<T> &col);
	NRmatrix<T> &setdiag(int jdiag, const NRvector<T> &diag);
	NRvector<T> sum(int dim) const;
	NRmatrix<T> transpose() const;
	NRmatrix<T> rotatecw() const;
	NRmatrix<T> rotateccw() const;
	NRmatrix<T> fliplr() const;
	NRmatrix<T> flipud() const;
	NRmatrix<T> hcat(const NRmatrix<T> &b);
	NRmatrix<T> vcat(const NRmatrix<T> &b);
	NRmatrix<T> repmat(int rows, int cols) const;
	NRmatrixview<T> gatherrows(const NRvector<int> locs);
	const NRmatrixview<T> gatherrows(const NRvector<int> locs) const;
	NRmatrix<T> gathercols(const NRvector<int> locs);
	NRmatrix<int> find() const; // return locs of nonzero entries for use with scatter/gather
	NRvector<T> gather(const NRmatrix<int> &locs) const;
	template <class U> NRmatrix<T> &scatter(const NRmatrix<int> &locs, const NRvector<U> &vals);
	template <class U> NRmatrix<T> &scatter(const NRmatrix<int> &locs, const U val);

	// slices all need const versions for const NRmatrices, nonconst for l-values
	NRmatrixview<T> slice(int rbeg, int rend, int cbeg, int cend);
	const NRmatrixview<T> slice(int rbeg, int rend, int cbeg, int cend) const;
	NRmatrixview<T> slice(char *star, int col);
	const NRmatrixview<T> slice(char *star, int col) const;
	NRmatrixview<T> slice(int row, char *star);
	const NRmatrixview<T> slice(int row, char *star) const;
	NRmatrixview<T> slice(char *star, int cbeg, int cend);
	const NRmatrixview<T> slice(char *star, int cbeg, int cend) const;
	NRmatrixview<T> slice(int rbeg, int rend, char *star);
	const NRmatrixview<T> slice(int rbeg, int rend, char *star) const;
	~NRmatrix();
	template <class U> friend class NRvector;
	template <class U> friend class NRvectorview;
	template <class U> friend class NRmatrixview;
};

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL) {
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, char *sep, ...) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)  {
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
    va_list ap;
    va_start(ap, sep);
    for (i=0;i<n;i++) for (j=0;j<m;j++) v[i][j] = va_arg(ap, T);
    va_end(ap);
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix<T> &rhs) : nn(rhs.nrows()), mm(rhs.ncols()) {
	int i,j,nel=mm*nn;
	v = (nn>0 ? new T*[nn] : NULL);
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T> template <class U>
NRmatrix<T>::NRmatrix(const NRmatrix<U> &rhs) : nn(rhs.nrows()), mm(rhs.ncols()) {
	int i,j,nel=mm*nn;
	v = (nn>0 ? new T*[nn] : NULL);
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) {
		v[i][j] = T(rhs[i][j]); // explicit conversion
	}
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs) {
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nrows() || mm != rhs.ncols()) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nrows();
			mm=rhs.ncols();
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const T rhs) {
	int i,j;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs;
	return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	THROW("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	THROW("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const { return nn; }
template <class T>
inline int NRmatrix<T>::ncols() const { return mm; }
template <class T>
inline int NRmatrix<T>::rend() const { return nn-1; }
template <class T>
inline int NRmatrix<T>::cend() const { return mm-1; }

template <class T>
NRmatrix<T>& NRmatrix<T>::resize(int newn, int newm, bool preserve) {
	int i,j,nel,nmin=MIN(nn,newn),mmin=MIN(mm,newm);
	if (newn != nn || newm != mm) {
		if (preserve) {
			T **vsave = v;
			v = newn>0 ? new T*[newn] : NULL;
			nel = newm*newn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< newn; i++) v[i] = v[i-1] + newm;
			for (i=0;i<nmin;i++) for (j=0;j<mmin;j++) v[i][j] = vsave[i][j];
			for (i=nmin;i<newn;i++) for (j=0;j<mmin;j++) v[i][j] = T(0);
			for (i=0;i<nmin;i++) for (j=mmin;j<newm;j++) v[i][j] = T(0);
			for (i=nmin;i<newn;i++) for (j=mmin;j<newm;j++) v[i][j] = T(0);
			if (vsave != NULL) {
				delete[] (vsave[0]);
				delete[] (vsave);
			}
			nn = newn;
			mm = newm;
		} else {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn = newn;
			mm = newm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
	}
	return *this;
}

template <class T>
NRmatrix<T>& NRmatrix<T>::assign(int newn, int newm, const T& a) {
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
	return *this;
}

template <class T>
T NRmatrix<T>::minval(int &rloc, int &cloc) const {
	int i,j, m=nrows(), n=ncols();
	T val = numeric_limits<T>::max();
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (v[i][j] < val) {
		val = v[i][j];
		rloc = i;
		cloc = j;
	}
	return val;
}
template <class T>
inline T NRmatrix<T>::minval() const {
	int rloc, cloc;
	return minval(rloc,cloc);
}

template <class T>
T NRmatrix<T>::maxval(int &rloc, int &cloc) const {
	int i,j, m=nrows(), n=ncols();
	T val = numeric_limits<T>::lowest();
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (v[i][j] > val) {
		val = v[i][j];
		rloc = i;
		cloc = j;
	}
	return val;
}
template <class T>
inline T NRmatrix<T>::maxval() const {
	int rloc, cloc;
	return maxval(rloc,cloc);
}

template <class T>
const NRvector<T> NRmatrix<T>::row(int irow) const {
	if (irow < 0 || irow > nn-1) THROW("impossible row index in NRmatrix::row");
	return NRvector<T>(mm,&v[irow][0]); // copies row of data
}


template <class T>
NRvectorview<T> NRmatrix<T>::row(int irow) {
	if (irow < 0 || irow > nn-1) THROW("impossible row index in NRmatrix::row");
	NRvectorview<T> ans;
	ans.nn = mm;
	ans.v = &v[irow][0];
	return ans;
}

template <class T>
NRvector<T> NRmatrix<T>::col(int jcol) const {
	if (jcol < 0 || jcol > mm-1) THROW("impossible col index in NRmatrix::col");
	NRvector<T> ans(nn);
	for (int i=0;i<nn;i++) ans[i] = v[i][jcol];
	return ans;
}

template <class T>
NRvector<T> NRmatrix<T>::diag(int jd) const {
	int i, ii;
	if (jd >= mm || -jd >= nn) THROW("impossible diag index in NRmatrix::diag");
	ii = (jd > 0 ? MAX(0,MIN(mm-jd,nn)) : MAX(0,MIN(mm,nn+jd)));
	NRvector<T> ans(ii);
	if (jd > 0) for (i=0;i<ii;i++) ans[i] = v[i][i+jd];
	else for (i=0;i<ii;i++) ans[i] = v[i-jd][i];
	return ans;
}

template <class T>
NRmatrix<T>&  NRmatrix<T>::setrow(int irow, const NRvector<T> &row) {
	if (irow < 0 || irow > nn-1) THROW("impossible row index in setrow");
	if (row.size() != mm) THROW("size mismatch in setrow");
	for (int j=0;j<mm;j++) v[irow][j] = row[j]; 
	return *this;
}

template <class T>
NRmatrix<T>&  NRmatrix<T>::setcol(int jcol, const NRvector<T> &col) {
	if (jcol < 0 || jcol > mm-1) THROW("impossible col index in setcol");
	if (col.size() != nn) THROW("size mismatch in setcol");
	for (int i=0;i<nn;i++) v[i][jcol] = col[i];
	return *this;
}

template <class T>
NRmatrix<T>&  NRmatrix<T>::setdiag(int jd, const NRvector<T> &diag) {
	int i, ii;
	if (jd >= mm || -jd >= nn) THROW("impossible diag index in NRmatrix::setdiag");
	ii = (jd > 0 ? MAX(0,MIN(mm-jd,nn)) : MAX(0,MIN(mm,nn+jd)));
	if (ii != diag.size()) THROW("size mismatch in setdiag");
	if (jd > 0) for (i=0;i<ii;i++) v[i][i+jd] = diag[i];
	else for (i=0;i<ii;i++) v[i-jd][i] = diag[i];
	return *this;
}


template <class T>
NRmatrixview<T> NRmatrix<T>::gatherrows(const NRvector<int> locs){
	int i, nl=locs.size();
	NRmatrixview<T> ans;
	ans.nn = nl;
	ans.mm = mm;
	ans.v = (nl > 0 ? new T*[nl] : NULL);
	for (i=0;i<nl;i++) ans.v[i] = v[locs[i]];
	return ans;
}

template <class T>
const NRmatrixview<T> NRmatrix<T>::gatherrows(const NRvector<int> locs) const{
	int i, nl=locs.size();
	NRmatrixview<T> ans;
	ans.nn = nl;
	ans.mm = mm;
	ans.v = (nl > 0 ? new T*[nl] : NULL);
	for (i=0;i<nl;i++) ans.v[i] = v[locs[i]];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::gathercols(NRvector<int> locs) {
	int jl,i,jcol, ml=locs.size();
	NRmatrix<T> ans(nn,ml);
	for (jl=0;jl<ml;jl++) {
		jcol = locs[jl];
		if (jcol < 0 || jcol > mm-1) THROW("impossible col index in NRmatrix::getcols");
		for (i=0;i<nn;i++) ans[i][jl] = v[i][jcol];
	}
	return ans;
}

template <class T>
NRmatrix<int> NRmatrix<T>::find() const {
	int i,j,n;
	n = 0;
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) if (v[i][j]) ++n;
	NRmatrix<int> ans(n,2);
	n = 0;
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) if (v[i][j]) {
		ans[n][0] = i;
		ans[n++][1] = j;
	}
	return ans;
}
template <class T>
NRvector<T> NRmatrix<T>::gather(const NRmatrix<int> &locs) const {
	int i, n = locs.nrows();
	NRvector<T> ans(n);
	for (i=0;i<n;i++) ans[i] = (operator[](locs[i][0]))[locs[i][1]]; // so that range checking is done if enabled
	return ans;
}
template <class T> template <class U>
NRmatrix<T> & NRmatrix<T>::scatter(const NRmatrix<int> &locs, const NRvector<U> &vals) {
	int i, n = locs.nrows();
	if (vals.size() != n) THROW("arguments in set must be same length");
	for (i=0;i<n;i++) (operator[](locs[i][0]))[locs[i][1]] = T(vals[i]); // so that range checking is done if enabled
	return *this;
}
template <class T> template <class U>
NRmatrix<T> & NRmatrix<T>::scatter(const NRmatrix<int> &locs, const U vals){
	int i, n = locs.nrows();
	if (vals.size() != n) THROW("arguments in set must be same length");
	for (i=0;i<n;i++) (operator[](locs[i][0]))[locs[i][1]] = T(vals); // so that range checking is done if enabled
	return *this;
}

template <class T>
T NRmatrix<T>::sum() const {
	int i,j;
	T ans = T(0);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans += v[i][j];
	return ans;
}

template <class T>
NRvector<T> NRmatrix<T>::sum(int dim) const {
	int i,j;
	T tmp;
	if (dim == 1) {
		NRvector<T> ans(mm);
		for (j=0;j<mm;j++) {
			tmp = T(0);
			for (i=0;i<nn;i++) tmp += v[i][j];
			ans[j] = tmp;
		}
		return ans;
	} else if (dim == 2) {
		NRvector<T> ans(nn);
		for (i=0;i<nn;i++) {
			tmp = T(0);
			for (j=0;j<mm;j++) tmp += v[i][j];
			ans[i] = tmp;
		}
		return ans;
	} else THROW("bad value of dim in matrix sum");
}

template <class T>
NRvectorview<T> NRmatrix<T>::asvec() {
	NRvectorview<T> ans;
	ans.nn = mm*nn;
	ans.v = &v[0][0];
	return ans;
}


/* template <class T>
NRvector<T> NRmatrix<T>::asvec() const {
	int i,len = mm*nn;
	T *ptr = &v[0][0];
	NRvector<T> ans(len);
	for (i=0;i<len;i++,ptr++) ans[i] = *ptr;
	return ans;
}*/

template <class T>
NRmatrix<T> NRmatrix<T>::transpose() const {
	int i,j;
	NRmatrix<T> ans(mm,nn);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans[j][i] = v[i][j];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::rotatecw() const {
	int i,j;
	NRmatrix<T> ans(mm,nn);
	for (i=0;i<mm;i++) for (j=0;j<nn;j++) ans[i][j] = v[nn-j-1][i];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::rotateccw() const {
	int i,j;
	NRmatrix<T> ans(mm,nn);
	for (i=0;i<mm;i++) for (j=0;j<nn;j++) ans[i][j] = v[j][mm-i-1];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::fliplr() const {
	int i,j;
	NRmatrix<T> ans(nn,mm);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans[i][j] = v[i][mm-j-1];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::flipud() const {
	int i,j;
	NRmatrix<T> ans(nn,mm);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans[i][j] = v[nn-i-1][j];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::hcat(const NRmatrix<T> &b) {
	int i,j,mb=b.ncols();
	if (nn != b.nrows()) THROW("incompatible shapes in hcat");
	NRmatrix<T> ans(nn,mm+mb);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans[i][j] = v[i][j];
	for (i=0;i<nn;i++) for (j=0;j<mb;j++) ans[i][j+mm] = b[i][j];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::vcat(const NRmatrix<T> &b) {
	int i,j, nb=b.nrows();
	if (mm != b.ncols()) THROW("incompatible shapes in vcat");
	NRmatrix<T> ans(nn+nb,mm);
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) ans[i][j] = v[i][j];
	for (i=0;i<nb;i++) for (j=0;j<mm;j++) ans[i+nn][j] = b[i][j];
	return ans;
}

template <class T>
NRmatrix<T> NRmatrix<T>::repmat(int nrow, int ncol) const {
	int ii,jj,i,j;
	NRmatrix<T> ans(nn*nrow,mm*ncol);
	for (ii=0;ii<nrow;ii++) for (jj=0;jj<ncol;jj++) {
		for (i=0;i<nn;i++) for (j=0;j<mm;j++) {
			ans[ii*nn+i][jj*mm+j] = v[i][j];
		}
	}
	return ans;
}

template <class T>
NRmatrix<T>::~NRmatrix() {
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

// NRmatrixview returned by slice "is a" NRmatrix that preserves data

template <class T>
class NRmatrixview : public NRmatrix<T> {
	using NRmatrix<T>::nn;
	using NRmatrix<T>::mm;
	using NRmatrix<T>::v;
public:
	NRmatrixview();
	NRmatrixview(const NRmatrixview<T> &rhs);
	NRmatrixview<T> & operator=(const NRmatrixview<T> &rhs);
	NRmatrixview<T> & operator=(const NRmatrix<T> &rhs);
	NRmatrixview<T> & operator=(const T rhs);
	~NRmatrixview(); 
};

template <class T>
NRmatrixview<T>::NRmatrixview() {}

template <class T>
NRmatrixview<T>::NRmatrixview(const NRmatrixview<T> &rhs) { // copy row ptrs, but not data
	nn = rhs.nn;
	mm = rhs.mm;
	v = (nn > 0 ? new T*[nn] : NULL);
	for (int i=0;i<nn;i++) v[i] = rhs.v[i];
};

template <class T>
NRmatrixview<T> & NRmatrixview<T>::operator=(const NRmatrixview<T> &rhs) {
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) THROW("incompatible sizes in slice assignment");
		for (i=0;i<nn;i++) for (j=0;j<mm;j++) v[i][j]=rhs[i][j];
	}
	return *this;
}

template <class T>
NRmatrixview<T> & NRmatrixview<T>::operator=(const NRmatrix<T> &rhs) {
	int i,j;
	if (nn != rhs.nn || mm != rhs.mm) THROW("incompatible sizes in slice assignment");
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) v[i][j]=rhs[i][j];
	return *this;
}

template <class T>
NRmatrixview<T> & NRmatrixview<T>::operator=(const T rhs) {
	int i,j;
	for (i=0;i<nn;i++) for (j=0;j<mm;j++) v[i][j]=rhs;
	return *this;
}

template <class T>
NRmatrixview<T>::~NRmatrixview() {
	if (v != NULL) delete [] v; // but don't delete v[0] which points to the data!
	v = NULL;
}

template <class T>
NRmatrixview<T> NRmatrix<T>::slice(int rbeg, int rend, int cbeg, int cend){
	int i;
	NRmatrixview<T> ans;
	ans.nn = MAX(0, rend - rbeg + 1);
	ans.mm = MAX(0, cend - cbeg + 1);
	ans.v = (ans.nn > 0 ? new T*[ans.nn] : NULL);
	for (i=0;i<ans.nn;i++) ans.v[i] = v[i+rbeg] + cbeg;
	return ans;
}
template <class T>
const NRmatrixview<T> NRmatrix<T>::slice(int rbeg, int rend, int cbeg, int cend) const {
	int i;
	NRmatrixview<T> ans;
	ans.nn = MAX(0, rend - rbeg + 1);
	ans.mm = MAX(0, cend - cbeg + 1);
	ans.v = (ans.nn > 0 ? new T*[ans.nn] : NULL);
	for (i=0;i<ans.nn;i++) ans.v[i] = v[i+rbeg] + cbeg;
	return ans;
}
template <class T>
inline NRmatrixview<T> NRmatrix<T>::slice(char *star, int col) {
	return slice(0,nn-1,col,col);
}
template <class T>
inline const NRmatrixview<T> NRmatrix<T>::slice(char *star, int col) const {
	return slice(0,nn-1,col,col);
}
template <class T>
inline NRmatrixview<T> NRmatrix<T>::slice(int row, char *star) {
	return slice(row,row,0,mm-1);
}
template <class T>
inline const NRmatrixview<T> NRmatrix<T>::slice(int row, char *star) const {
	return slice(row,row,0,mm-1);
}
template <class T>
inline NRmatrixview<T> NRmatrix<T>::slice(char *star, int cbeg, int cend) {
	return slice(0,nn-1,cbeg,cend);
}
template <class T>
inline const NRmatrixview<T> NRmatrix<T>::slice(char *star, int cbeg, int cend) const {
	return slice(0,nn-1,cbeg,cend);
}
template <class T>
inline NRmatrixview<T> NRmatrix<T>::slice(int rbeg, int rend, char *star) {
	return slice(rbeg,rend,0,mm-1);
}
template <class T>
inline const NRmatrixview<T> NRmatrix<T>::slice(int rbeg, int rend, char *star) const {
	return slice(rbeg,rend,0,mm-1);
}

// end NRmatrix and related

template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n]) {
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d() {
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

// basic type names (redefine if your bit lengths don't match)

typedef int32_t Int; // 32 bit integer
typedef uint32_t Uint;

typedef int64_t Llong; // 64 bit integer
typedef uint64_t Ullong;

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;

typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;


// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
static turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

// sleep in milliseconds (supposed to be portable?) WHP 5/28/18
static void sleep(int ms) {
	chrono::milliseconds timespan(ms);
	this_thread::sleep_for(timespan);
}

template <class T>
NRvector<T> range(T start, T stop, T incr=T(1)) {
	T val = start;
	incr = (stop >= start ? abs(incr) : -abs(incr));
	int i,n = 1 + int((stop - start)/incr + T(0.5)); // nearest integer, including the endpoint
	NRvector<T> ans(n);
	for (i=0;i<n;i++) {ans[i] = val; val += incr;} // MUNG
	return ans;
}

template <class T>
NRvector<T> matmul(NRmatrix<T> &a, NRvector<T> &b) {
	// matrix vector multiplication
	int i,k,m=a.nrows(), p=a.ncols();
	T sum;
	if (p != b.size()) THROW("incompatible sizes in mat*vec matmul");
	NRvector<T> ans(m);
	for (i=0;i<m;i++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[i][k]*b[k];
		ans[i] = sum;
	}
	return ans;
}

template <class T>
NRvector<T> matmul(NRmatrix<T> &a, char *transpose, NRvector<T> &b) {
	// matrix^T vector multiplication
	int i,k,m=a.ncols(), p=a.nrows();
	T sum;
	if (p != b.size()) THROW("incompatible sizes in mat^T * vec matmul");
	NRvector<T> ans(m);
	for (i=0;i<m;i++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[k][i]*b[k];
		ans[i] = sum;
	}
	return ans;
}

template <class T>
NRvector<T> matmul(NRvector<T> &a, NRmatrix<T> &b) {
	// vector matrix multiplication
	int j,k, n=b.ncols(), p=b.nrows();
	T sum;
	if (p != a.size()) THROW("incompatible sizes in vec*mat matmul");
	NRvector<T> ans(n);
	for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[k]*b[k][j];
		ans[j] = sum;
	}
	return ans;
}

template <class T>
NRvector<T> matmul(NRvector<T> &a, NRmatrix<T> &b, char *transpose) {
	// vector matrix^T multiplication
	int j,k, n=b.nrows(), p=b.ncols();
	T sum;
	if (p != a.size()) THROW("incompatible sizes in vec * mat^T matmul");
	NRvector<T> ans(n);
	for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[k]*b[j][k];
		ans[j] = sum;
	}
	return ans;
}

template <class T>
NRmatrix<T> matmul(NRmatrix<T> &a, NRmatrix<T> &b) {
	// matrix multiplication
	int i,j,k,m=a.nrows(), n=b.ncols(), p=a.ncols();
	T sum;
	if (p != b.nrows()) THROW("incompatible sizes in matmul 1");
	NRmatrix<T> ans(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[i][k]*b[k][j];
		ans[i][j] = sum;
	}
	return ans;
}

template <class T>
NRmatrix<T> matmul(NRmatrix<T> &a, char *transpose, NRmatrix<T> &b) {
	// matrix multiplication a^T * b
	int i,j,k,m=a.ncols(), n=b.ncols(), p=a.nrows();
	T sum;
	if (p != b.nrows()) THROW("incompatible sizes in matmul 2");
	NRmatrix<T> ans(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[k][i]*b[k][j];
		ans[i][j] = sum;
	}
	return ans;
}

template <class T>
NRmatrix<T> matmul(NRmatrix<T> &a, NRmatrix<T> &b, char *transpose) {
	// matrix multiplication a * b^T
	int i,j,k,m=a.nrows(), n=b.nrows(), p=a.ncols();
	T sum;
	if (p != b.ncols()) THROW("incompatible sizes in matmul 3");
	NRmatrix<T> ans(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[i][k]*b[j][k];
		ans[i][j] = sum;
	}
	return ans;
}

template <class T>
NRmatrix<T> matmul(NRmatrix<T> &a, char *t1, NRmatrix<T> &b, char *t2) {
	// matrix multiplication
	int i,j,k,m=a.ncols(), n=b.nrows(), p=a.nrows();
	T sum;
	if (p != b.ncols()) THROW("incompatible sizes in matmul 4");
	NRmatrix<T> ans(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) {
		sum = T(0);
		for (k=0;k<p;k++) sum += a[k][i]*b[j][k];
		ans[i][j] = sum;
	}
	return ans;
}

template <class T>
T dotproduct(const NRvector<T> &a, const NRvector<T> &b) {
	int i, m=a.size();
	if (m != b.size()) THROW("incompatible sizes in dotproduct");
	T sum = T(0);
	for (i=0;i<m;i++) sum += a[i]*b[i];
	return sum;
}

template <class T>
NRmatrix<T> outerproduct(const NRvector<T> &a, const NRvector<T> &b) {
	int i,j, m=a.size(), n=b.size();
	NRmatrix<T> ans(m,n);
	for (i=0;i<m;i++) for (j=0;j<n;j++) ans[i][j] = a[i]*b[j];
	return ans;
}

// element-wise binary ops
#define NewVecMatBinOp(OPEQ,OPEQNAME,OPNAME) \
template <class T> \
NRvector<T> & OPEQNAME (NRvector<T> &a, const NRvector<T> &b) { \
	Int i, m = a.size(); \
	if (m != b.size()) THROW("incompatible sizes in vector op"); \
	for (i=0; i<m; i++) a[i] OPEQ b[i]; \
	return a; \
} \
template <class T> \
NRvector<T> & OPEQNAME (NRvector<T> &a, T s) { \
	Int i, m = a.size(); \
	for (i=0; i<m; i++) a[i] OPEQ s; \
	return a; \
} \
template <class T> \
NRvector<T> OPNAME (const NRvector<T> &a, const NRvector<T> &b) \
{ return NRvector<T>(a) OPEQ b; } \
template <class T> \
NRvector<T> OPNAME (const NRvector<T> &a, T s) \
{ return NRvector<T>(a) OPEQ s; } \
template <class T> \
NRmatrix<T> & OPEQNAME (NRmatrix<T> &a, const NRmatrix<T> &b) { \
	Int i,j, m = a.nrows(), n = a.ncols(); \
	if (m != b.nrows() || n != b.ncols()) THROW("incompatible sizes in matrix op"); \
	for (i=0;i<m;i++) for (j=0;j<n;j++) a[i][j] OPEQ b[i][j]; \
	return a; \
} \
template <class T> \
NRmatrix<T> & OPEQNAME (NRmatrix<T> &a, T s) { \
	Int i,j, m = a.nrows(), n = a.ncols(); \
	for (i=0;i<m;i++) for (j=0;j<n;j++) a[i][j] OPEQ s; \
	return a; \
} \
template <class T> \
NRmatrix<T> OPNAME (const NRmatrix<T> &a, const NRmatrix<T> &b) \
{ return NRmatrix<T>(a) OPEQ b; } \
template <class T> \
NRmatrix<T> OPNAME (const NRmatrix<T> &a, T s) \
{ return NRmatrix<T>(a) OPEQ s; }

NewVecMatBinOp(+=,operator+=,operator+)
NewVecMatBinOp(-=,operator-=,operator-)
NewVecMatBinOp(*=,operator*=,operator*)
NewVecMatBinOp(/=,operator/=,operator/)
NewVecMatBinOp(%=,operator%=,operator%)
NewVecMatBinOp(&=,operator&=,operator&)
NewVecMatBinOp(|=,operator|=,operator|)
NewVecMatBinOp(^=,operator^=,operator^)

// element-wise unary ops that allow two args
#define NewVecMatUnaryOp2(NAME,OP) \
template <class T> \
NRvector<T>& NAME(const char *inplace, NRvector<T> &a) { \
	int i, n=a.size(); \
	for (i=0;i<n;i++) a[i] = OP(a[i]); \
	return a; \
} \
template <class T> \
NRvector<T> NAME(const NRvector<T> &a) { \
	NRvector<T> ans(a); \
	return NAME("inplace", ans); \
} \
template <class T> \
NRmatrix<T>& NAME(const char *inplace, NRmatrix<T> &a) { \
	int i,j, n=a.nrows(), m=a.ncols(); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) a[i][j] = OP(a[i][j]); \
	return a; \
} \
template <class T> \
NRmatrix<T> NAME(const NRmatrix<T> &a) { \
	NRmatrix<T> ans(a); \
	return NAME("inplace", ans); \
}

NewVecMatUnaryOp2(sqr,SQR); // deprecated: use SQR directly
NewVecMatUnaryOp2(operator-,-);

// element-wise unary ops that don't allow two args
#define NewVecMatUnaryOp(NAME,OP) \
template <class T> \
NRvector<T> NAME(const NRvector<T> &a) { \
	int i, n=a.size(); \
	NRvector<T> ans(n); \
	for (i=0;i<n;i++) ans[i] = OP(a[i]); \
	return ans; \
} \
template <class T> \
NRmatrix<T> NAME(const NRmatrix<T> &a) { \
	int i,j, n=a.nrows(), m=a.ncols(); \
	NRmatrix<T> ans(n,m); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) ans[i][j] = OP(a[i][j]); \
	return ans; \
} 

NewVecMatUnaryOp(operator!,!);
NewVecMatUnaryOp(operator~,~);


// logical binary ops
#define NewVecMatBinaryLogicalOp(NAME,OP) \
template <class T> \
NRvector<unsigned char> NAME(const NRvector<T> &a, const NRvector<T> &b) { \
	Int i, m = a.size(); \
	if (m != b.size()) THROW("incompatible sizes in vector logical op"); \
	NRvector<unsigned char> ans(m); \
	for (i=0;i<m;i++) ans[i] = (a[i] OP b[i]); \
	return ans; \
} \
template <class T> \
NRmatrix<unsigned char> NAME(const NRmatrix<T> &a, const NRmatrix<T> &b) { \
	Int i,j, n = a.nrows(), m = a.ncols(); \
	if (n != b.nrows() || m != b.ncols()) THROW("incompatible sizes in matrix logical op"); \
	NRmatrix<unsigned char> ans(n,m); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) ans[i][j] = (a[i][j] OP b[i][j]); \
	return ans; \
} \
template <class T> \
NRvector<unsigned char> NAME(const NRvector<T> &a, T b) { \
	Int i, m = a.size(); \
	NRvector<unsigned char> ans(m); \
	for (i=0;i<m;i++) ans[i] = (a[i] OP b); \
	return ans; \
} \
template <class T> \
NRmatrix<unsigned char> NAME(const NRmatrix<T> &a, T b) { \
	Int i,j, n = a.nrows(), m = a.ncols(); \
	NRmatrix<unsigned char> ans(n,m); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) ans[i][j] = (a[i][j] OP b); \
	return ans; \
}

NewVecMatBinaryLogicalOp(operator==,==)
NewVecMatBinaryLogicalOp(operator!=,!=)
NewVecMatBinaryLogicalOp(operator<,<)
NewVecMatBinaryLogicalOp(operator<=,<=)
NewVecMatBinaryLogicalOp(operator>,>)
NewVecMatBinaryLogicalOp(operator>=,>=)
NewVecMatBinaryLogicalOp(operator&&,&&)
NewVecMatBinaryLogicalOp(operator||,||)

// threaded common functions

#define NewThreadedFunction(FUNC) \
template <class T> NRvector<T>& FUNC(char *inplace, NRvector<T> &a) { \
	int i,n=a.size(); \
	for (i=0;i<n;i++) a[i] = FUNC(a[i]); \
	return a; \
} \
template <class T> NRvector<T> FUNC(const NRvector<T> &a) { \
	NRvector<T> ans(a); \
	FUNC("inplace",ans); \
	return ans; \
} \
template <class T> NRmatrix<T>& FUNC(char *inplace, NRmatrix<T> &a) { \
	int i,j,n=a.nrows(),m=a.ncols(); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) a[i][j] = FUNC(a[i][j]); \
	return a; \
} \
template <class T> NRmatrix<T> FUNC(const NRmatrix<T> &a) { \
	NRmatrix<T> ans(a); \
	FUNC("inplace",ans); \
	return ans; } 

NewThreadedFunction(sin)
NewThreadedFunction(cos)
NewThreadedFunction(tan)
NewThreadedFunction(asin)
NewThreadedFunction(acos)
NewThreadedFunction(atan)
NewThreadedFunction(cosh)
NewThreadedFunction(sinh)
NewThreadedFunction(tanh)
NewThreadedFunction(exp)
NewThreadedFunction(log)
NewThreadedFunction(log10)
NewThreadedFunction(sqrt)
NewThreadedFunction(abs)
NewThreadedFunction(fabs)
NewThreadedFunction(ceil)
NewThreadedFunction(floor)

#define NewTwoArgThreadedFn(OP) \
template<class T> \
NRvector<T> OP(const NRvector<T> &a, const NRvector<T> &b) { \
	Int i,n=a.size(); \
	if (b.size() != n) THROW("incompatible sizes in NewTwoArgThreadedFn"); \
	NRvector<T> ans(a); \
	for (i=0;i<n;i++) ans[i] = OP(a[i],b[i]); \
	return ans; \
} \
template<class T> \
NRmatrix<T> OP(const NRmatrix<T> &a, const NRmatrix<T> &b) { \
	Int i,j,n=a.nrows(),m=a.ncols(); \
	if (b.nrows() != n || b.ncols() != m) THROW("incompatible sizes in SIGN"); \
	NRmatrix<T> ans(a); \
	for (i=0;i<n;i++) for (j=0;j<m;j++) ans[i][j] = OP(a[i][j],b[i][j]); \
	return ans; \
}
NewTwoArgThreadedFn(SIGN)
NewTwoArgThreadedFn(MAX)
NewTwoArgThreadedFn(MIN)

// following must be static to keep the linker from seeing them as duplicated
static char nr_out_float_format[] = "%.5f ";
static char nr_out_fixed_format[] = "%d ";
static char nr_out_ullong_format[] = "%016llx ";
static char nr_out_unsignedfixed_format[] = "%x ";
static char nr_in_float_format[] = "%lf";
static char nr_in_fixed_format[] ="%d";
static char *nr_out_format(Doub x) {return nr_out_float_format;}
static char *nr_out_format(Int x) {return nr_out_fixed_format;}
static char *nr_out_format(Uint x) {return nr_out_unsignedfixed_format;}
static char *nr_out_format(Ullong x) {return nr_out_ullong_format;}
static char *nr_out_format(Uchar x) {return nr_out_fixed_format;}
static char *nr_in_format(Doub x) {return nr_in_float_format;}
static char *nr_in_format(Int x) {return nr_in_fixed_format;}
static char *nr_in_format(Uchar x) {return nr_in_fixed_format;}

template <class T>
void printmat(const NRmatrix<T> &a, const char *format=NULL) {
	int i,j,m=a.nrows(), n=a.ncols();
	if (format == NULL) format = nr_out_format(T(1));
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) printf(format,a[i][j]);
		printf("\n");
	}
	printf("\n");
}

template <class T>
void printmat(const char *precomment, const NRmatrix<T> &a, const char *format=NULL) {
	int i,j,m=a.nrows(), n=a.ncols();
	if (format == NULL) format = nr_out_format(T(1));
	printf("%s (%d by %d matrix)\n",precomment,m,n);
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) printf(format,a[i][j]);
		printf("\n");
	}
	printf("\n");
}

template <class T>
void printvec(const NRvector<T> &a, const char *format=NULL) {
	int j, n=a.size();
	if (format == NULL) format = nr_out_format(T(1));
	for (j=0;j<n;j++) printf(format,a[j]);
	printf("\n");
}

template <class T>
void printvec(const char *precomment, const NRvector<T> &a, const char *format=NULL) {
	int j, n=a.size();
	if (format == NULL) format = nr_out_format(T(1));
	printf("%s (length %d vector)\n",precomment, n);
	for (j=0;j<n;j++) printf(format,a[j]);
	printf("\n");
}

template <class T>
void loadmat(NRmatrix<T> &X, char *filename, char* format=NULL) {
	Int i,j,n=X.nrows(), m=X.ncols();
	if (format == NULL) format = nr_in_format(T(1));
	FILE *INP = fopen(filename,"rb");
	if (!INP) THROW("bad file in loadmat");
	for (i=0;i<n;i++) for (j=0;j<m;j++)
		if (fscanf(INP,format,&X[i][j]) != 1) THROW("bad read in loadmat");
	fclose(INP);
}

template <class T>
void loadvec(NRvector<T> &y, char *filename, char* format=NULL) {
	Int i,n=y.size();
	if (format == NULL) format = nr_in_format(T(1));
	FILE *INP = fopen(filename,"rb");
	if (!INP) THROW("bad file in loadvec");
	for (i=0;i<n;i++)
		if (fscanf(INP,format,&y[i]) != 1) THROW("bad read in loadvec");
	fclose(INP);
}

template <class T>
void dumpmat(NRmatrix<T> &X, char *filename, char* format=NULL) {
	Int i,j,n=X.nrows(), m=X.ncols();
	if (format == NULL) format = nr_out_format(T(1));
	FILE *OUTP = fopen(filename,"wb");
	if (!OUTP) THROW("bad file in dumpmat");
	for (i=0;i<n;i++) {
		for (j=0;j<m;j++) fprintf(OUTP,format,X[i][j]);
		fprintf(OUTP,"\n");
	}
	fclose(OUTP);
}

template <class T>
void dumpvec(NRvector<T> &y, char *filename, char* format=NULL) {
	Int i,n=y.size();
	if (format == NULL) format = nr_out_format(T(1));
	FILE *OUTP = fopen(filename,"wb");
	if (!OUTP) THROW("bad file in dumpvec");
	for (i=0;i<n;i++) fprintf(OUTP,format,y[i]);
	fclose(OUTP);
}

//#endif /*__INTELLISENSE__*/
#endif /* _NR3_H_ */
