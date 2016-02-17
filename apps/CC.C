// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"
#include <cmath>

/*
alg: 
- assign start times
- while (exists unvisited vertices) :
    - add in any centers that should start on this step
    - expand frontier using CAS to acquire
- now, have a map from id -> cluster_id
- run parallel union-find. only care about inter-component edges
- 
*/


/**
  Generates a integer sequence of shifts, where shifts[i] is the start
  of the vertices to take on round i.
 **/
uintT *generateShifts(intT n, double beta) {
  //only need (ln n)/beta levels
  uintT maxLevel = min<uintT>(n+1,2+ceil(log(n)/beta));
  cout << "maxLevel is : " << maxLevel << endl;
  uintT *betas = newA(uintT, maxLevel);
  {parallel_for(intT i=0; i<maxLevel; i++) {
    betas[i] = floor(exp(i*beta));
  }}
  uintT last = sequence::plusScan(betas, betas, maxLevel);
  return betas;
}


struct CC_F {
  uintE* components;
  CC_F(uintE* _components) : components(_components) {}
  inline bool update (uintE s, uintE d) { //Update
    if (components[d] == UINT_E_MAX) {
      components[d] = components[s];
      return 1; 
    }
    else return 0;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic version of Update
    return (CAS(&components[d],UINT_E_MAX,components[s]));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return (components[d] == UINT_E_MAX); }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long n = GA.n;

  double beta = 0.1;
  uintT *shifts = generateShifts(n, beta);

  // note: can randomly permute here---leaving out for now

  //creates components array, initialized to maxint
  uintE* components = newA(uintE,n);
  parallel_for(long i=0;i<n;i++) components[i] = UINT_E_MAX;

  intT* flags = newA(intT, n+1);

  vertexSubset Frontier(n); 
  long totalVisited = 0;
  intT round = 0;
  while (!Frontier.isEmpty() || totalVisited < n) {

    intT start = shifts[round];
    intT end = min<uintT>(shifts[round+1],n); 
    intT size = end - start;

    parallel_for(intT i=0;i<size;i++) {
      flags[i] = (components[start+i] == UINT_E_MAX);
    }
    flags[size] = 0;
    intT numAdded = sequence::plusScan(flags,flags,size+1);

    // no easy way to meld two frontiers right now. One option is to 
    // pre-allocate an array and then repeatedly write into it (creating a
    // vertexSubset out of it on each round), but it might get turned into a
    // dense version at some point, which prevents us from using it as in the 
    // original CC algorithm
    // For now, I will just create a new frontier explicitly and fix this with
    // a more optimized version soon. 
    
    if (!Frontier.isEmpty()) {
      Frontier.toSparse();
    }
    long prevFrontierSize = Frontier.numNonzeros();
    long nextFrontierSize = prevFrontierSize + numAdded;
    uintE* nextF = newA(uintE, nextFrontierSize); 
    // copy old frontier
    parallel_for(intT i = 0; i < prevFrontierSize; i++) {
      nextF[i] = Frontier.s[i];
    }

    parallel_for(intT i = 0; i < size; i++) {
      intT offset = flags[i];
      if(flags[i+1] != offset) {
        nextF[prevFrontierSize+offset] = start+i;
        components[start+i] = start+i;
      }
    }

    Frontier.del();
    Frontier = vertexSubset(n, nextFrontierSize, nextF);

    totalVisited += Frontier.numNonzeros();

    vertexSubset output = edgeMap(GA, Frontier, CC_F(components),GA.m/20);
    Frontier.del();
    Frontier = output;
    round += 1;
  }

  long ic_edges = 0;
  parallel_for (int i = 0; i < n; i++) {
    long our_ic = 0;
    for (int j = 0; j < GA.V[i].getOutDegree(); j++) {
      if (components[i] != components[GA.V[i].getOutNeighbor(j)]) {
        our_ic++;
      }
    }
    writeAdd<long>(&ic_edges, our_ic);
  }

  for (int i = 0; i < n; i++) {
    cout << "components[i] = " << components[i] << endl;
  }

  cout << "m = " << GA.m << endl;
  cout << "ic edges = " << ic_edges << endl;

  cout << "done" << endl;
  
  Frontier.del();
  free(components); 
}
