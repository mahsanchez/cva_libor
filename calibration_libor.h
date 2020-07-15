#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>

#include "mkl_lapacke.h"

#include "curvepoint.h"
#include "leastsquares.h"

//#define DEBUG_CAPSTRIPPING 10

std::vector<double> tenor = {
    //3M, 6M, 9M, 1Y, 1Y3M, 1Y6M, 1Y9M, 2Y, 2Y3M, 2Y6M, 2Y9M, 3Y, 3Y3M, 3Y6M, 3Y9M, 4Y, 4Y3M, 4Y6M, 4Y9M, 5Y
    0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.50, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
};

std::vector<double> yearFractions = {
     0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25278, 0.25278, 0.26111,
     0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278
};

std::vector<double> bonds = {
    0.9947527, 0.9892651, 0.9834984, 0.9774658, 0.9712884, 0.9648035, 0.9580084, 0.9509789, 0.9440868, 0.9369436, 0.9295484, 0.9219838,
    0.9145031, 0.9068886, 0.8990590, 0.8911017, 0.8833709, 0.8754579, 0.8673616
};

// MArket Cap Volatility
std::vector<double> capvols = {
    0.1641, 0.1641, 0.1641, 0.1765, 0.1889, 0.2013, 0.2137, 0.2162, 0.2186, 0.2211, 0.2235, 0.2223, 0.2212, 0.2200, 0.2188, 0.2173, 0.2158, 0.2142, 0.2127
};

// Fitted Volatility
std::vector<double> vols {
    0.1641, 0.1641, 0.1641, 0.2015, 0.2189, 0.2365, 0.2550, 0.2212, 0.2255, 0.2298, 0.2341, 0.2097, 0.2083, 0.2077, 0.2051, 0.2007, 0.1982, 0.1959, 0.1938
};

std::vector<std::vector<double>> rho(19, std::vector<double>(19));

// TODO Add Least Square Solution

#ifndef __INTERVAL_BISECTION_H
#define __INTERVAL_BISECTION_H

template<typename T>
double interval_bisection(T f, const std::random_device::result_type entropy)
{
    std::mt19937 gen(entropy);
    static const auto lower_bound = 0.01;
    static const auto upper_bound =  0.80;
    std::uniform_real_distribution<> dis(lower_bound, upper_bound);

    auto pos_pt = dis(gen);
    auto neg_pt = dis(gen);

    while (f(pos_pt) < 0.0)
        pos_pt = dis(gen);

    while (f(neg_pt) > 0.0)
        neg_pt = dis(gen);

    static const auto about_zero_mag = 1E-6;
    for (;;)
    {
        const auto mid_pt = (pos_pt + neg_pt)/2.0;
        const auto f_mid_pt = f(mid_pt);
        if (fabs(f_mid_pt)  < about_zero_mag)
            return mid_pt;

        if (f_mid_pt >= 0.0)
            pos_pt = mid_pt;
        else
            neg_pt = mid_pt;
    }
}

#endif


class cpl {
public:
       cpl(double _bond, double _yearFraction, double _strike, double _x, double _rate) :
         bond(_bond), yearFraction(_yearFraction), strike(_strike), x(_x), rate(_rate)
       {}

       inline double d(double vol) {
            double result = std::log(rate/x) + 0.5*vol*vol * yearFraction;
            result /= vol * std::sqrt(yearFraction);
            return result;
       }

       double operator()(double vol) {
           double d1 = d(vol);
           double d2 = d1 - vol * std::sqrt(yearFraction);
           double result = rate * N( d1 ) - strike * N( d2 );
           result *= bond;
           return result;
       }

       double N(double val) {
           double result = 0;
           vdCdfNorm( 1, &val, &result );
           return result;
       }

private:
    double bond;
    double x;
    double yearFraction;
    double strike;
    double rate;
};


class CplVolCalibration {
public:
    CplVolCalibration(std::vector<double> &yearFraction_, std::vector<double> &bonds_, std::vector<double> &capvols_, double expiry_) :
     yearFraction(yearFraction_), bonds(bonds_), capvols(capvols_), expiry(expiry_), dtau(0.25) {
        strikes.resize(bonds.size(), 0.0);
        cplvols.resize(bonds.size(), 0.0);
        rates.resize(bonds.size(), 0.0);
    }

    // At the Money Cap Strike [S􏴩(T0, T3m, Ti)]
    void caps_strikes() {
        // calculate ATM strikes for caps
        std::vector<double> terms(strikes.size());

        // yearFraction * DF
        terms[0] = 0;
        for (int i = 1; i < terms.size(); i++) {
            terms[i] = yearFraction[i] * bonds[i];
        }

        // cumulative sum of terms
        std::partial_sum(terms.begin() + 1, terms.end(), strikes.begin() + 1, std::plus<double>());

        // strikes calculation
        strikes[0] = 0.0;
        for (int i = 1; i < strikes.size(); i++) {
            strikes[i] = (bonds[0] - bonds[i]) / strikes[i];
        }
    }

    void caplet_volatility_stripping() {
        //
        int size = expiry/dtau;

        // Interpolation for 3M, 6M, 9M from 1Y store those values in capvols[]
        cplvols[0] = 0.1641;     // 3M
        cplvols[1] = 0.1641;     // 6M
        //cplvols[2] = 0.1641;     // 9M
        //cplvols[3] = 0.1641;     // 1Y

        // Calculate Rates and Interpolation
        for (int i = 1; i < bonds.size(); i++) {
            rates[i] = ( bonds[i-1]/bonds[i] -1 ) / yearFraction[i];
        }

        //Price the caps 1Y3M (Onwards) up to 5Y
        std::vector<double> caps(size);
        std::vector<double> caplet_prices(size);

        for (int i = 2; i < size; i++) {
            double result = 0.0;
            for (int j = 1; j < i; j++) {
                result += cpl(bonds[j], yearFraction[j], strikes[j], strikes[j], rates[j])(capvols[j]);
            }
            caps[i-1] = result;
        }

        std::random_device rd;
        static const auto entropy = rd();

        //Caculate Implied Caplet Volatilities
        for (int i = 1; i < size-1; i++) {
            // caplet_price
            double caplet_price = 0.0;
            for (int j = 1; j < i; j++) {
                caplet_price += cpl(bonds[j], yearFraction[j], strikes[j], strikes[j], rates[j])(cplvols[j]);
            }

            // Caplet Pricer to be used for Implied Volatility with Root Finding
            cpl cpl_pricer(bonds[i], yearFraction[i], strikes[i], strikes[i], rates[i]);

            auto cpl_functor = [&] (double vol) {
#ifdef DEBUG_CAPSTRIPPING
                std::cout << i << " cpl_functor " << caps[i] - caplet_price << " " << caps[i] << " " << caplet_price << " " << vol  << std::endl;
#endif
                return caps[i] - caplet_price - cpl_pricer(vol);
            };

            // root finding using interval bisection algorithm
            cplvols[i] = interval_bisection(cpl_functor, entropy);

#ifdef DEBUG_CAPSTRIPPING
            std::cout << i << " capvol vs cplvol" << capvols[i] << " " << cplvols[i]  << std::endl;
#endif
        }

    }

    /* Volatility Fitting - solve optimization problem  argmin sum[cplvols - fitvols]^2 using least square methods */
    void volatility_fitting() {
        int size = expiry/dtau;

        /*
        // instvol = phi * [ (a + b*(T - t)) * exp(-c* (T - t)) + d] * 1 {t < T}
        for (int i = 0; i < size; i++) {
            fitvols[i] = phi * ( (a + b*(yearFraction[i] - t)) * std::exp(-c* (yearFraction[i]  - t)) + d );
        }
        //Find a,b, c, d using Least Square methods (intel MKL dgels)
        //int info = LAPACKE_dgels( LAPACK_ROW_MAJOR, 'N', m, n, nrhs, &cplvols[0], lda, &fitvols[0], ldb );
        */

        const double phi = 1.0;
        double a = -0.1;
        double b = 0.5;
        double c = 1.0;
        double d = 0.1;

        curveFitting();

        // instvols is a 1-d vector instvol = phi * [ (a + b*(T - t)) * exp(-c* (T - t)) + d] * 1 {t < T}
        instvols.resize(size, std::vector<double>(size));

        // convert from instanteneous volatility to LMM volatilities using piece wise 2-d vector
        for (int i = 1; i < size; i++) {
            for (int j = 0; j < i; j++) {
                double term = (tenor[i-1] - tenor[j]);
                instvols[i][j] = phi * ( (a + b*term) * std::exp(-c*term) + d );
            }
        }

    }

    /*
     * QuantLib Least Squares Method implementation
    // http://mikejuniperhill.blogspot.com/2016/04/quantlib-least-squares-method.html
     */
    void curveFitting() {
        // 2nd degree polynomial least squares fitting for a curve
        // create observed market rates for a curve
        Array independentValues(20);
        independentValues[0] = 0.1641; independentValues[1] = 0.1641; independentValues[2] = 0.1641;
        independentValues[3] = 0.2015; independentValues[4] = 0.2189; independentValues[5] = 0.2365;
        independentValues[6] = 0.2550; independentValues[7] = 0.2212; independentValues[8] = 0.2255;
        independentValues[9] = 0.2298; independentValues[10] = 0.2341; independentValues[11] = 0.2097;
        independentValues[12] = 0.2083; independentValues[13] = 0.2077; independentValues[14] = 0.2051;
        independentValues[15] = 0.2007; independentValues[16] = 0.1982; independentValues[17] = 0.1959;
        independentValues[18] = 0.1938; independentValues[19] = 0.1925;
        //
        const double phi = 1.0;
        double a = 0.1;
        double b = 0.1;
        double c = 0.1;
        double d = 0.1;

        auto calibration_functor = [&](double T, double t) {
            return phi * ( (a + b*(T - t)) * std::exp(-c*(T - t)) + d );
        };

        // create corresponding curve points to be approximated
        std::vector<boost::shared_ptr<CurvePoint>> curvePoints;
        for (int i = 0; i < independentValues.size(); i++) {
            curvePoints.push_back(boost::shared_ptr<CurvePoint>(
                  new CurvePoint(calibration_functor(tenor[i], 0.0) )));
        }

        //
        // create container for function pointers for calculating rate approximations
        std::vector<boost::function<Real(const Array&)>> dependentFunctions;
        //
        // for each curve point object, bind function pointer to operator() and add it into container
        for (unsigned int i = 0; i < curvePoints.size(); i++)
        {
            dependentFunctions.push_back(boost::bind(&CurvePoint::operator(), curvePoints[i], _1));
        }
        // perform least squares fitting and print optimized coefficients
        LeastSquares leastSquares(boost::shared_ptr<OptimizationMethod>(new LevenbergMarquardt), 10000, 1000, 1E-09, 1E-09, 1E-09);
        NoConstraint parametersConstraint;
        Array initialParameters(4, 0.1);
        Array coefficients = leastSquares.Fit(independentValues, dependentFunctions, initialParameters, parametersConstraint);

        a = initialParameters[0];
        b = initialParameters[1];
        c = initialParameters[2];
        d = initialParameters[3];

        for (int i = 0; i < independentValues.size(); i++) {
            std::cout << phi * ( (a + b*(tenor[i] - 0.0)) * std::exp(-c*(tenor[i] - 0.0)) + d ) << std::endl;
        }

        for (unsigned int i = 0; i < coefficients.size(); i++)
        {
            std::cout << coefficients[i] << std::endl;
        }
    }


    void correlation() {
        double beta = 0.1;
        int size = expiry/dtau;
        rho.resize(size, std::vector<double>(size));

        for (int i = 1; i < size; i++) {
            for (int j = 0; j < size; j++) {
                rho[i-1][j] = std::exp( -beta * (tenor[i] - tenor[j]) );
            }
        }
    }

    void calibrate() {
        caps_strikes();
        caplet_volatility_stripping();
        volatility_fitting();
        correlation();
    }

private:
    double dtau;
    double expiry;
    std::vector<double> &yearFraction;
    std::vector<double> &bonds;
    std::vector<double> &capvols;
    std::vector<double> rates;
    std::vector<double> strikes;
    std::vector<double> cplvols;
    std::vector<std::vector<double>> instvols;
    std::vector<std::vector<double>> rho;

};
