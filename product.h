#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <memory>
#include <cmath>
#include <iostream>

#include "q_numerics.h"

#include <ql/quantlib.hpp>
using namespace QuantLib;
using namespace std;

#define DEBUG 0

// Exposure Points
std::vector<double> exposure_timepoints = {
        1,2,3,4,5, 15, 22, 29, 30, 36, 50, 64, 76, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 780, 900, 1020, 1140, 1260, 1380, 1500, 1720, 1840, 2140, 2520, 2880, 3240, 3600, 4800, 7200
};


/*
 * Pricing Instrument Interest Rate Swap IRS
 */

/*
 * VainillaInterestRateSwap
 * Valuation = Sum[0..T] (S - K) *  𝛼j * DF(t0, tj)
 *
 * 􏰷 𝛼j is the day-count fraction over which the coupon accrues.
􏰷 * B(t0 , tj ) is the value (at t0 ) of a discount factor maturing at time tj .
􏰷 * K is the fixed rate of the swap.
 * 􏰷S is the par swap rate where
 *
􏰷 The fixed and floating leg frequencies and day count bases are assumed to be the same for simplicity.
 */
struct InterestRateSwap {
    InterestRateSwap(std::vector<double> &pricing_points_, std::vector<double> &floating_schedule_,  std::vector<double> &fixed_schedule_, double notional_, double K_, double expiry_, double dtau_) :
            pricing_points(pricing_points_), floating_schedule(floating_schedule_), fixed_schedule(fixed_schedule_), notional(notional_), K(K_), expiry(expiry_), dtau(dtau_)
    {}

    std::vector<double> &pricing_points;
    std::vector<double> &floating_schedule;
    std::vector<double> &fixed_schedule;
    double notional;
    double K;
    double dtau;
    double expiry;
};

/**
 * Interpolated ZCB Curve
 */
class InterpolatedFICurve {
public :
    InterpolatedFICurve(std::vector<double>& dates, std::vector<double>& values) :
            dates_(dates), values_(values),  interpolation_ (LinearInterpolation(dates.begin(), dates.end(), values.begin())) {}

    double operator() (double t) {
        return  interpolation_(t, true);
    }

private:
    std::vector<double> dates_;
    std::vector<double> values_;
    LinearInterpolation interpolation_;
};

/**
 * Spreads interpolated Curve
 */
class InterpolatedSpreadsCurve {
public :
    InterpolatedSpreadsCurve(std::vector<double>& dates, std::vector<double>& values, double bps_) :
            dates_(dates), values_(values),  interpolation_ (LinearInterpolation(dates.begin(), dates.end(), values.begin())), bps(bps_) {}

    double operator() (double t) {
        return interpolation_(t, true) * bps;
    }
private:
    double bps;
    std::vector<double> dates_;
    std::vector<double> values_;
    LinearInterpolation interpolation_;
};


/*
 * Survival Probability TermStructure implementation
 * CDS Bootstrapping JPMorgan Methodology
 * VB code at http://mikejuniperhill.blogspot.com/2014/08/bootstrapping-default-probabilities.html
*/

class SurvivalProbCurve {
public:
    SurvivalProbCurve(std::vector<double>& timepoints_, InterpolatedSpreadsCurve &spreadsCurve_, InterpolatedFICurve &zcbCurve_, double recovery) :
            timepoints(timepoints_), spreads(spreadsCurve_),  zcb(zcbCurve_), interpolation_(LinearInterpolation(timepoints.begin(), timepoints.end(), timepoints.begin()))
    {
        probabilities.resize(timepoints.size());
        bootstrap(probabilities, timepoints, recovery);
        interpolation_ =  LinearInterpolation(timepoints.begin(), timepoints.end(), probabilities.begin());
    }

    double operator() (double timepoint) {
        return interpolation_(timepoint);
    }

    void bootstrap(std::vector<double>& p, std::vector<double>& timepoints, double recovery) {
        double loss = 1 - recovery;
        double term, terms, divider, term1, term2;

        for (int i = 0; i < timepoints.size(); i++) {
            if (i == 0) {
                p[0] = 1.0;
            }
            else if (i == 1) {
                p[1] = loss / (spreads(timepoints[1])* (timepoints[1] - timepoints[0]) + loss);
            }
            else {
                terms = 0.0;
                double spread = spreads(timepoints[i]);
                for (int j = 1; j < i; j++) {
                    double dtau = timepoints[j] - timepoints[j-1];
                    term = loss * p[j-1];
                    term -= (loss + dtau * spread) * p[j];
                    term *= zcb(timepoints[j]);
                    terms += term;
                }

                double dtau = timepoints[i] - timepoints[i-1];
                divider = loss + dtau * spread;
                divider *= zcb(timepoints[i]);
                term1 = terms/divider;
                term2 = p[i-1] * loss;
                term2 /= (loss + dtau * spread);
                p[i] = term1 + term2;
            }
        }
    }

private:
    std::vector<double> probabilities;
    std::vector<double>& timepoints;
    InterpolatedFICurve &zcb;
    InterpolatedSpreadsCurve &spreads;
    LinearInterpolation interpolation_;
};


