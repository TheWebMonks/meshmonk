#ifndef NEIGHBOURFINDER_HPP
#define NEIGHBOURFINDER_HPP

#include "meshmonk/global.hpp"
#include "meshmonk/profiling.hpp"
#include <Eigen/Dense>
#include <nanoflann.hpp>

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
    MatDynInt; // matrix MxN of type unsigned int
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;
typedef Eigen::Matrix<float, Eigen::Dynamic, registration::NUM_FEATURES>
    FeatureMat;
typedef Eigen::Matrix<float, Eigen::Dynamic, 3> Vec3Mat;

namespace registration {

// Process-wide default for the NeighbourFinder warm-start cache
// (ADR-007). Each instance copies this value on construction, so
// callers must set it BEFORE invoking register functions. Tests flip
// it to false to verify cached and uncached paths give identical
// output.
extern bool g_neighbour_caching_default;

template <typename VecMatType> class NeighbourFinder {
  /*
  GOAL
  This class searches for the k nearest neighbours in 'inSourcePoints' for
  each element in the 'inQueriedPoints' set. It outputs the indices of each
  neighbour and the squared (!) distances between each element of
  'inQueriedPoints' and its neighbours.

  INPUT
  -inQueriedPoints:
  -inSourcePoints:

  PARAMETERS
  -numNeighbours(= 3): number of nearest neighbours
  -leafSize(= 15): should be between 5 and 50 or so

  OUTPUT
  -outNeighbourIndices
  -outNeighbourSquaredDistances.
  */

public:
  // NeighbourFinder();
  ~NeighbourFinder(); // destructor

  void set_source_points(const VecMatType *const inQueriedPoints);
  void set_queried_points(const VecMatType *const _inSourcePoints);
  MatDynInt get_indices() const { return _outNeighbourIndices; }
  MatDynFloat get_distances() const { return _outNeighbourSquaredDistances; }
  void set_parameters(const size_t numNeighbours);
  void update();

  // Disable the cross-iteration warm-start cache. Used in tests to verify
  // cached and uncached code paths produce identical results.
  void set_caching_enabled(bool enabled) { _cachingEnabled = enabled; }

protected:
private:
  // # Inputs
  const VecMatType *_inQueriedPoints = NULL;
  const VecMatType *_inSourcePoints = NULL;

  // # Outputs
  MatDynInt _outNeighbourIndices;
  MatDynFloat _outNeighbourSquaredDistances;
  // # User parameters

  // # Internal Data structures
  nanoflann::KDTreeEigenMatrixAdaptor<VecMatType> *_kdTree = NULL;

  // # Warm-start cache across update() calls (ADR-007 D1)
  MatDynInt _cachedNeighbourIndices;
  const VecMatType *_cachedSourcePtr = NULL;
  size_t _cachedNumSourceElements = 0;
  bool _cachingEnabled = g_neighbour_caching_default;

  // # Interal parameters
  size_t _numDimensions = 0;
  size_t _numSourceElements = 0;
  size_t _numQueriedElements = 0;
  size_t _numNeighbours = 3;
  size_t _leafSize = 15;
};

/*
See
http://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file
for an explanation why we instantiate our templated class with every matrix type
we will use here.
*/
//
template class NeighbourFinder<FeatureMat>;
template class NeighbourFinder<Vec3Mat>;
//

template <typename VecMatType> NeighbourFinder<VecMatType>::~NeighbourFinder() {
  // destructor
  if (_kdTree != NULL) {
    delete _kdTree;
    _kdTree = NULL;
  }
}

template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_source_points(
    const VecMatType *const inSourcePoints) {
  // # Cache invalidation (ADR-007 D3): same pointer + same row count ⇒
  // retain cache; cached indices still refer to the same mesh vertices
  // whose positions may have changed in place. Any other change drops
  // the cache to avoid dereferencing stale indices.
  const size_t newNumSourceElements =
      static_cast<size_t>(inSourcePoints->rows());
  if (inSourcePoints != _cachedSourcePtr ||
      newNumSourceElements != _cachedNumSourceElements) {
    _cachedNeighbourIndices.resize(0, 0);
  }
  _cachedSourcePtr = inSourcePoints;
  _cachedNumSourceElements = newNumSourceElements;

  // # Set input
  _inSourcePoints = inSourcePoints;

  // # Update internal parameters
  _numDimensions = _inSourcePoints->cols();
  _numSourceElements = _inSourcePoints->rows();

  // # Update internal data structures
  // ## The kd-tree has to be rebuilt.
  {
#ifdef MESHMONK_PROFILING
    auto _t_build = g_profiler.scoped("NeighbourFinder::build_kdtree");
#endif
    if (_kdTree != NULL) {
      delete _kdTree;
      _kdTree = NULL;
    }
    _kdTree = new nanoflann::KDTreeEigenMatrixAdaptor<VecMatType>(
        *_inSourcePoints, _leafSize);
    _kdTree->index->buildIndex();
  }
}

template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_queried_points(
    const VecMatType *const inQueriedPoints) {
  // # Set input
  _inQueriedPoints = inQueriedPoints;

  // # Update internal parameters
  _numQueriedElements = _inQueriedPoints->rows();

  // # Adjust internal data structures
  // ## The indices and distance matrices have to be resized.
  _outNeighbourIndices.setZero(_numQueriedElements, _numNeighbours);
  _outNeighbourSquaredDistances.setZero(_numQueriedElements, _numNeighbours);
}

template <typename VecMatType>
void NeighbourFinder<VecMatType>::set_parameters(const size_t numNeighbours) {
  // # Check if what user requests, changes the parameter value
  bool parameterChanged = false;
  if (_numNeighbours != numNeighbours) {
    parameterChanged = true;
  }

  // # Set parameter
  _numNeighbours = numNeighbours;

  // # Resize the output matrices if the parameter is changed
  if (parameterChanged == true) {
    _outNeighbourIndices.setZero(_numQueriedElements, _numNeighbours);
    _outNeighbourSquaredDistances.setZero(_numQueriedElements, _numNeighbours);
    // Cache shape is numQueried × old-k; invalidate on k change.
    _cachedNeighbourIndices.resize(0, 0);
  }
}

template <typename VecMatType> void NeighbourFinder<VecMatType>::update() {
#ifdef MESHMONK_PROFILING
  auto _t = g_profiler.scoped("NeighbourFinder::update");
#endif

  // # Warm-start cache (ADR-007): cache is usable iff it was written by
  // a prior update() with the same shape. Fresh-search path on mismatch.
  const bool useCache =
      _cachingEnabled &&
      static_cast<size_t>(_cachedNeighbourIndices.rows()) == _numQueriedElements &&
      static_cast<size_t>(_cachedNeighbourIndices.cols()) == _numNeighbours;

  // ADR-007 D2 — must match the eps passed to SearchParams below.
  const float searchEps = 0.0001f;

  // # Query the kd-tree
  // ## Loop over the queried features
  // ### Hoist buffers above loop (D2): first allocation happens here;
  // ### resize() inside the loop is a no-op after iteration 0.
  unsigned int i = 0;
  unsigned int j = 0;
  std::vector<float> queriedFeature;
  std::vector<size_t> neighbourIndices;
  std::vector<float> neighbourSquaredDistances;
  nanoflann::KNNResultSet<float> knnResultSet(_numNeighbours);

  // ### Execute loop
  for (; i < _numQueriedElements; ++i) {
    // resize() is a no-op after iteration 0 (vectors stay at target capacity)
    queriedFeature.resize(_numDimensions);
    neighbourIndices.resize(_numNeighbours);
    neighbourSquaredDistances.resize(_numNeighbours);
    knnResultSet = nanoflann::KNNResultSet<float>(_numNeighbours);

    // ### Initiliaze the knnResultSet
    knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);

    // ### convert input features to 'queriedFeature' std::vector structure
    // ### (required by nanoflann's kd-tree).
    for (j = 0; j < _numDimensions; ++j) {
      queriedFeature[j] = (*_inQueriedPoints)(i, j);
    }

    // ### Warm-start: compute the max distance from the query to the
    // previous iteration's neighbour indices at CURRENT positions. This
    // is a valid upper bound on the true k-th distance (ADR-007 D1).
    // Overwrite the worstDist slot set by init() — safe because
    // KNNResultSet::addPoint only reads dists[0..count-1] while
    // count < capacity; dists[capacity-1] is the worstDist sentinel
    // until count reaches capacity.
    if (useCache) {
      float ceiling = 0.0f;
      for (size_t c = 0; c < _numNeighbours; ++c) {
        const int cachedIdx = _cachedNeighbourIndices(i, c);
        float d = 0.0f;
        for (size_t dim = 0; dim < _numDimensions; ++dim) {
          const float diff = queriedFeature[dim] -
                             (*_inSourcePoints)(cachedIdx, dim);
          d += diff * diff;
        }
        if (d > ceiling) {
          ceiling = d;
        }
      }
      neighbourSquaredDistances[_numNeighbours - 1] =
          ceiling * (1.0f + searchEps);
    }

    // tree_query — kd-tree nearest-neighbour search (dominant cost bucket)
    {
#ifdef MESHMONK_PROFILING
      auto _t_query = g_profiler.scoped("NeighbourFinder::tree_query");
#endif
      _kdTree->index->findNeighbors(
          knnResultSet, &queriedFeature[0],
          nanoflann::SearchParams(32, searchEps, true));
    }

    for (j = 0; j < _numNeighbours; ++j) {
      _outNeighbourIndices(i, j) = neighbourIndices[j];
      _outNeighbourSquaredDistances(i, j) = neighbourSquaredDistances[j];
    }
  }

  if (_cachingEnabled) {
    _cachedNeighbourIndices = _outNeighbourIndices;
  }
} // end k_nearest_neighbours()

} // namespace registration

#endif // NEIGHBOURFINDER_HPP
