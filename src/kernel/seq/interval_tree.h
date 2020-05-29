/* Copyright 2019 The Spin-Scenario Authors. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#pragma once
#include <vector>
#include <algorithm>

namespace ssl {
namespace seq {

template<class T, typename K = std::size_t>
class Interval {
 public:
  K start;
  K stop;
  T value;
  Interval(K s, K e, const T &v)
      : start(s),
        stop(e),
        value(v) {
  }
};

template<class T, typename K>
K intervalStart(const Interval<T, K> &i) {
  return i.start;
}

template<class T, typename K>
K intervalStop(const Interval<T, K> &i) {
  return i.stop;
}

template<class T, typename K>
std::string &operator<<(std::string &out, Interval<T, K> &i) {
  out << "interval(" << i.start << ", " << i.stop << "): " << i.value;
  return out;
}

template<class T, typename K = std::size_t>
class IntervalStartSorter {
 public:
  bool operator()(const Interval<T, K> &a, const Interval<T, K> &b) {
    return a.start < b.start;
  }
};

template<class T, typename K = std::size_t>
class IntervalTree {

 public:
  typedef Interval<T, K> interval;
  typedef std::vector<interval> intervalVector;
  typedef IntervalTree<T, K> intervalTree;

  intervalVector intervals;
  intervalTree *left;
  intervalTree *right;
  K center;

  IntervalTree<T, K>(void)
      : left(NULL),
        right(NULL),
        center(0) {
  }

  IntervalTree<T, K>(const intervalTree &other)
      : left(NULL),
        right(NULL) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = new intervalTree(*other.left);
    }
    if (other.right) {
      right = new intervalTree(*other.right);
    }
  }

  IntervalTree<T, K> &operator=(const intervalTree &other) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = new intervalTree(*other.left);
    } else {
      if (left)
        delete left;
      left = NULL;
    }
    if (other.right) {
      right = new intervalTree(*other.right);
    } else {
      if (right)
        delete right;
      right = NULL;
    }
    return *this;
  }

  IntervalTree<T, K>(intervalVector &ivals, std::size_t depth = 16,
                     std::size_t minbucket = 64, K leftextent = 0,
                     K rightextent = 0, std::size_t maxbucket = 512)
      : left(NULL),
        right(NULL) {

    --depth;
    IntervalStartSorter<T, K> intervalStartSorter;
    if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
      std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
      intervals = ivals;
    } else {
      if (leftextent == 0 && rightextent == 0) {
        // sort intervals by start
        std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
      }

      K leftp = 0;
      K rightp = 0;
      K centerp = 0;

      if (leftextent || rightextent) {
        leftp = leftextent;
        rightp = rightextent;
      } else {
        leftp = ivals.front().start;
        std::vector<K> stops;
        stops.resize(ivals.size());
        transform(ivals.begin(), ivals.end(), stops.begin(),
                  intervalStop<T, K>);
        rightp = *max_element(stops.begin(), stops.end());
      }

      //centerp = ( leftp + rightp ) / 2;
      centerp = ivals.at(ivals.size() / 2).start;
      center = centerp;

      intervalVector lefts;
      intervalVector rights;

      for (typename intervalVector::iterator i = ivals.begin();
           i != ivals.end(); ++i) {
        interval &interval = *i;
        if (interval.stop < center) {
          lefts.push_back(interval);
        } else if (interval.start > center) {
          rights.push_back(interval);
        } else {
          intervals.push_back(interval);
        }
      }

      if (!lefts.empty()) {
        left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
      }
      if (!rights.empty()) {
        right = new intervalTree(rights, depth, minbucket, centerp, rightp);
      }
    }
  }

  void findOverlapping(K start, K stop, intervalVector &overlapping) const {
    if (!intervals.empty() && !(stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin();
           i != intervals.end(); ++i) {
        const interval &interval = *i;
        if (interval.stop > start && interval.start < stop) {
          overlapping.push_back(interval);
        }
      }
    }

    if (left && start <= center) {
      left->findOverlapping(start, stop, overlapping);
    }

    if (right && stop >= center) {
      right->findOverlapping(start, stop, overlapping);
    }

  }

  void findContained(K start, K stop, intervalVector &contained) const {
    if (!intervals.empty() && !(stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin();
           i != intervals.end(); ++i) {
        const interval &interval = *i;
        if (interval.start >= start && interval.stop <= stop) {
          contained.push_back(interval);
        }
      }
    }

    if (left && start <= center) {
      left->findContained(start, stop, contained);
    }

    if (right && stop >= center) {
      right->findContained(start, stop, contained);
    }

  }

  ~IntervalTree(void) {
    // traverse the left and right
    // delete them all the way down
    if (left) {
      delete left;
    }
    if (right) {
      delete right;
    }
  }

};
} /* namespace seq */
} /* namespace ssl */
