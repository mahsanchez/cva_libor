#include <utility>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>

#include "mkl.h"
#include "mkl_vsl.h"


/**
 * percentile
 */

double percentile(std::vector<double> &vectorIn, double percent)
{
    auto nth = vectorIn.begin() + (percent*vectorIn.size())/100;
    std::nth_element(vectorIn.begin(), nth, vectorIn.end());
    return *nth;
}

/*
 * Prefix Sum
 */
void prefix_sum(double *prefixSum, double *arr, int n)
{
    prefixSum[0] = arr[0];
    // Adding present element with previous element
    for (int i = 1; i < n; i++) {
        prefixSum[i] = prefixSum[i - 1] + arr[i];
    }
}


class linear_interpolator {
public:
    linear_interpolator() {
    }

    linear_interpolator(const std::vector<double>& x, const std::vector<double>& values) {
        initialize(x, values);
    }

    void initialize(const std::vector<double>& x, const std::vector<double>& values)
    {
        for (int i = 0; i < values.size(); i++) {
            points.push_back( {x[i], values[i]} );
        }
    }

    double operator() (double x) const {
         return find(x);
    }

    void add(double x, double value) {
        points.push_back( {x, value} );
    }

    std::vector<std::pair<double, double>>& getPoints() {
        std::vector<std::pair<double, double>> &result = points;
        return points;
    }

    double find(double x) const {
        if (points.size() == 0) {
            return 0.0;
        }

        //the element already exists no interpolation is needed
        auto pointIter = std::find_if(points.begin(), points.end(), [x](std::pair<double, double> point ) {
            if (point.first == x) {
                return true;
            }
            return false;
        });

        if (pointIter != std::end(points))  {
           return pointIter->second;
        }

        //Define a lambda that returns true if the x value of a point pair is < the caller's x value
        auto lessThan =
                [](const std::pair<double, double>& point, double x)
                {return point.first < x;};

        //Find the first table entry whose value is >= caller's x value
        auto iter =
                std::lower_bound(points.cbegin(), points.cend(), x, lessThan);

        //If the caller's X value is greater than the largest
        //X value in the table, we can't interpolate.
        if(iter == points.cend()) {
            return (points.cend() - 1)->second;
        }

        //If the caller's X value is less than the smallest X value in the table,
        //we can't interpolate.
        if(iter == points.cbegin() and x <= points.cbegin()->first) {
            return points.cbegin()->second;
        }

        //We can interpolate!
        double upperX{iter->first};
        double upperY{iter->second};
        double lowerX{(iter - 1)->first};
        double lowerY{(iter - 1)->second};
        double deltaY{upperY - lowerY};
        double deltaX{upperX - lowerX};

        return lowerY + ((x - lowerX)/ deltaX) * deltaY;
    }

private:
    std::vector<std::pair<double, double>> points;
};



/*
 * matrix transposition rectangular matrix
Inspired by the Wikipedia - Following the cycles algorithm description, I came up with following C++ implementation:
 https://stackoverflow.com/questions/9227747/in-place-transposition-of-a-matrix
 https://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose/?ref=rp

#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <algorithm> // std::swap (until C++11)
#include <vector>

template<class RandomIterator>
void transpose(RandomIterator first, RandomIterator last, int m)
{
    const int mn1 = (last - first - 1);
    const int n   = (last - first) / m;
    std::vector<bool> visited(last - first);
    RandomIterator cycle = first;
    while (++cycle != last) {
        if (visited[cycle - first])
            continue;
        int a = cycle - first;
        do  {
            a = a == mn1 ? mn1 : (n * a) % mn1;
            std::swap(*(first + a), *cycle);
            visited[a] = true;
        } while ((first + a) != cycle);
    }
}

int test()
{
    int a[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    transpose(a, a + 8, 4);
    std::copy(a, a + 8, std::ostream_iterator<int>(std::cout, " "));
}

**/