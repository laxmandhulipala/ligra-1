#ifndef VERTEXSUBSET_H
#define VERTEXSUBSET_H

struct vertexSubset {
  long n, m;
  uintE* s;
  bool* d;
  bool isDense;

  // make a singleton vertex in range of n
vertexSubset(long _n, intE v) 
: n(_n), m(1), d(NULL), isDense(0) {
  s = newA(uintE,1);
  s[0] = v;
}
  
  //empty vertex set
vertexSubset(long _n) : n(_n), m(0), d(NULL), s(NULL), isDense(0) {}
  // make vertexSubset from array of vertex indices
  // n is range, and m is size of array
vertexSubset(long _n, long _m, uintE* indices) 
: n(_n), m(_m), s(indices), d(NULL), isDense(0) {}
  // make vertexSubset from boolean array, where n is range
vertexSubset(long _n, bool* bits) 
: n(_n), d(bits), s(NULL), isDense(1)  {
  m = sequence::sum(bits,_n); }
  // make vertexSubset from boolean array giving number of true values
vertexSubset(long _n, long _m, bool* bits) 
: n(_n), m(_m), s(NULL), d(bits), isDense(1)  {}

  // delete the contents
  void del(){
    if (d != NULL) free(d);
    if (s != NULL) free(s);
  }
  long numRows() { return n; }
  long numNonzeros() { return m; }
  bool isEmpty() { return m==0; }

  // converts to dense but keeps sparse representation if there
  void toDense() {
    if (d == NULL) {
      d = newA(bool,n);
      {parallel_for(long i=0;i<n;i++) d[i] = 0;}
      {parallel_for(long i=0;i<m;i++) d[s[i]] = 1;}
    }
    isDense = true;
  }

  // converts to sparse but keeps dense representation if there
  void toSparse() {
    if (s == NULL) {
      _seq<uintE> R = sequence::packIndex<uintE>(d,n);
      if (m != R.n) {
	cout << "bad stored value of m" << endl; 
	abort();
      }
      s = R.A;
    }
    isDense = false;
  }
  // check for equality
  bool eq (vertexSubset& b) {
    toDense();
    b.toDense();
    bool* c = newA(bool,n);
    {parallel_for (long i=0; i<b.n; i++) 
	c[i] = (d[i] != b.d[i]);}
    bool equal = (sequence::sum(c,n) == 0);
    free(c);
    return equal;
  }

  void print() {
    if (isDense) {
      cout << "D:";
      for (long i=0;i<n;i++) if (d[i]) cout << i << " ";
      cout << endl;
    } else {
      cout << "S:";
      for (long i=0; i<m; i++)	cout << s[i] << " ";
      cout << endl;
    }
  }
};

//*****VERTEX FUNCTIONS*****

//Note: this is the optimized version of vertexMap which does not
//perform a filter
template <class F>
void vertexMap(vertexSubset V, F add) {
  long n = V.numRows(), m = V.numNonzeros();
  if(V.isDense) {
    {parallel_for(long i=0;i<n;i++)
	if(V.d[i]) add(i);}
  } else {
    {parallel_for(long i=0;i<m;i++)
	add(V.s[i]);}
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options {DENSE, DENSE_FORWARD};

#endif
