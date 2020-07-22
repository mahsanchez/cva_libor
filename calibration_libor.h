#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>

#include "mkl_lapacke.h"

#include "curvepoint.h"
#include "leastsquares.h"

#include <ql/math/solvers1d/bisection.hpp>

#define DEBUG_CAPSTRIPPING 10
#define DEBUG_IMPLIEDCAPLETVOL 10
#define DEBUG_FITTEDVOL 10

// number of days since (Ti - T0)/360
std::vector<double> tenor = {
    //3M, 6M, 9M, 1Y, 1Y3M, 1Y6M, 1Y9M, 2Y, 2Y3M, 2Y6M, 2Y9M, 3Y, 3Y3M, 3Y6M, 3Y9M, 4Y, 4Y3M, 4Y6M, 4Y9M, 5Y
    0.2514, 0.5028, 0.7583, 1.0139, 1.2639, 1.5167, 1.7722, 2.0278, 2.2778, 2.5306, 2.7861, 3.0417, 3.2944, 3.5472, 3.8083, 4.0611, 4.3139, 4.5667, 4.8194, 5.0722
};

std::vector<double> yearFractions = {
     0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25278, 0.25278, 0.26111,
     0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278
};

std::vector<double> bonds = { //zcb
    0.9947527, 0.9892651, 0.9834984, 0.9774658, 0.9712884, 0.9648035, 0.9580084, 0.9509789, 0.9440868, 0.9369436, 0.9295484, 0.9219838,
    0.9145031, 0.9068886, 0.8990590, 0.8911017, 0.8833709, 0.8754579, 0.8673616, 0.8581725
};

// MArket Cap Volatility
std::vector<double> capvols = {
    0.0, 0.1641, 0.1641, 0.1641, 0.1765, 0.1889, 0.2013, 0.2137, 0.2162, 0.2186, 0.2211, 0.2235, 0.2223, 0.2212, 0.2200, 0.2188, 0.2173, 0.2158, 0.2142, 0.2127
};

/*
 * Price a CAPLET using Black Formula
 */
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

/*
 * Price a CAPLET using Black Formula
 */
double caplet_price(double zcb, double yearFraction, double strike, double x, double rate, double vol) {
    cpl caplet(zcb, yearFraction, strike, x, rate);
    return caplet(vol);
}

/**
 * Implied Caplet Volatility Problem LMM in Practice Algorithm 7.1
 */
double impliedVolProb(double zcb, double yearFraction, double strike, double x, double fwd_rate, double vol, double cap, double cplet_price) {
    return (cap - cplet_price - caplet_price(zcb, yearFraction,  strike, x, fwd_rate, vol));
}

/**
 * Stripped Caplet Volatility
 */

class StrippedCapletVolatility {
public:
    StrippedCapletVolatility(std::vector<double> &yearFraction_, std::vector<double> &zcb_, std::vector<double> &capvols_, double expiry_, double dtau_) :
    yearFraction(yearFraction_), zcb(zcb_), capvols_mkt(capvols_), expiry(expiry_), dtau(dtau_) {
        strikes.resize(zcb.size(), 0.0);
    }

    // Calculate Forward Swap Rate (ATM strikes) [S􏴩(T0, T3m, Ti)] LMM in Practice chapter 7 p73 Table 7.4
    void caps_strikes() {
        // calculate ATM strikes for caps
        std::vector<double> terms(strikes.size());

        // yearFraction * DF
        terms[0] = 0;
        for (int i = 1; i < terms.size(); i++) {
            terms[i] = yearFraction[i] * zcb[i];
        }

        // cumulative sum of terms
        std::partial_sum(terms.begin(), terms.end(), strikes.begin(), std::plus<double>());

        // strikes calculation
        strikes[0] = 0.0;
        for (int i = 1; i < strikes.size(); i++) {
            strikes[i] = (zcb[0] - zcb[i]) / strikes[i];
        }

#ifdef DEBUG_CAPSTRIPPING
        std::cout << "Forward Swap Rate (ATM strikes) [S(T0, T3m, Ti)]" << std::endl;
        for (int i = 0; i < strikes.size(); i++) {
            std::cout << tenor[i]  << " " << strikes[i] << std::endl;
        }
#endif
    }

    // Implied Caplet Volatility - LMM in Practice Algorithm 7.1
    void caplet_volatility_stripping(std::vector<double> cplvols, std::vector<double> &spot_rates) {
        int size = zcb.size();

        // calculate spot_rates
        spot_rates.resize(zcb.size());

        // Calculate Spot Rates (LIBOR Rate covering) from a ZCB Curve, LMM in Practice p70 equation 7.2
        for (int i = 1; i < zcb.size(); i++) {
            spot_rates[i] = ( (zcb[i-1]/zcb[i]) - 1 ) * (1.0 / yearFraction[i-1]);
        }

#ifdef DEBUG_IMPLIEDCAPLETVOL
        std::cout << "Spot Rates" << std::endl;
        for (int i = 0; i < spot_rates.size(); i++) {
            std::cout << tenor[i]  << " " << spot_rates[i] << std::endl;
        }
#endif

        //Price the caps 1Y3M (Onwards) up to 5Y
        std::vector<double> caps_quotes(size);
        std::vector<double> caplet_prices(size);

        std::vector<double> acc(size);
        for (int i = 1; i < size; i++) {
            acc[i] = caplet_price(zcb[i], yearFraction[i], strikes[i], strikes[i], spot_rates[i], capvols_mkt[i]);
        }
        std::partial_sum(acc.begin(), acc.end(), caps_quotes.begin());

#ifdef DEBUG_IMPLIEDCAPLETVOL
        std::cout << "CAP Quotes Prices" << std::endl;
        for (int i = 0; i < caps_quotes.size(); i++) {
            std::cout << tenor[i]  << " " << caps_quotes[i] << std::endl;
        }
#endif
        // 1D Root Finding Solver using Bisection Method
        Bisection root_solver;

        // TODO - Interpolation for 3M, 6M, 9M from 1Y store those values in capvols[]
        // LMM in Practice chapter 7 p76 Algorithm 7.1
        cplvols.resize(size);
        cplvols[1] = 0.1641;     // 6M
        //Caculate Implied Caplet Volatilities
        for (int i = 2; i < cplvols.size(); i++) {
            double cplet_price = 0.0; // caplet_price
            for (int j = 1; j < i; j++) {
                cplet_price += caplet_price(zcb[j], yearFraction[j], strikes[j], strikes[j], spot_rates[j],cplvols[j]);
            }
            //calculate implied caplet volatility
            boost :: function < Real ( Volatility )> impliedVolFunc ;   // setup a boost function
            // bind the boost function to all market parameters , keep vol as variant
            impliedVolFunc = boost :: bind (& impliedVolProb , zcb[i], yearFraction[i], strikes[i], strikes[i], spot_rates[i], _1, caps_quotes[i], cplet_price);
            // root finding using interval bisection algorithm
            cplvols[i] = root_solver.solve(impliedVolFunc, 1E-6, 0.1600, 0.1500, 0.25);
        }

#ifdef DEBUG_IMPLIEDCAPLETVOL
        // LMM in Practice chapter 7 p77 Table 7.5
        std::cout << "Caplet Volatilites Stripped from cap volatilities" << std::endl;
        for (int i = 0; i < cplvols.size(); i++) {
            std::cout << tenor[i]  << " " << capvols[i] << " " << cplvols[i] << std::endl;
        }
#endif
    }

private:
    double dtau;
    double expiry;
    std::vector<double> &yearFraction;
    std::vector<double> &zcb;
    std::vector<double> &capvols_mkt;
    std::vector<double> strikes;
    std::vector<double> capletvols = {
            0.1641, 0.1641, 0.1641, 0.2015, 0.2189, 0.2365, 0.2550, 0.2212, 0.2255, 0.2298,
            0.2341, 0.2097, 0.2083, 0.2077, 0.2051, 0.2007, 0.1982, 0.1959, 0.1938
    };
};

/**
 *
 *  CALIBRATION
 */

/* Parametric Calibration LMM in Practice Chapter 9 p 158 Equation (9.12) */

inline double parametric_calibration(double t, double a, double b, double c, double d, double Ti) {
    double val = a + (b + c*(Ti - t)) * std::exp(-d*(Ti - t));
    return val;
};

inline double squared_parametric_calibration(double t, double a, double b, double c, double d, double Ti) {
    double val = parametric_calibration(t, a, b, c, d, Ti);
    return val * val;
};

/*
 * Optimization Contraints
 * Constraint the a,b,c,d domain such that parametric function a + (b + c*(Ti - t)) * std::exp(-d*(Ti - t)) shows an exponential decae
 */
class FittingVolConstraint : public Constraint {
private:
    class Impl : public Constraint::Impl {
    public:
        bool test(const Array& params) const {
            if ( (params[0] + params[1]) <= 0) return false;
            if (params[3] <= 0.0) return false;
            return true;
        }
        Array upperBound(const Array& params) const {
            return Array(params.size(),2.0);
        }
        Array lowerBound(const Array& params) const {
            Array lower_bound(params.size(), 0.0);
            lower_bound[1] = -0.2;
            lower_bound[3] = 0.9;
            return lower_bound;
        }
    };
public:
    FittingVolConstraint()
            : Constraint(ext::shared_ptr<Constraint::Impl>(new FittingVolConstraint::Impl)) {}
};


/* Caplet Volatility & Volatility Fitting */

class CplVolCalibration {
public:
    CplVolCalibration(std::vector<double> &yearFraction_, std::vector<double> &bonds_, std::vector<double> &capletvols_, double expiry_) :
       yearFraction(yearFraction_), bonds(bonds_), cplvols(capletvols_), expiry(expiry_), dtau(0.25) {}

    /*
     * Volatility Fitting - solve optimization problem  argmin sum[cplvols - fitvols]^2 using least square methods
     * LMM in practice Chapter 9
     * QuantLib Least Squares Method implementation http://mikejuniperhill.blogspot.com/2016/04/quantlib-least-squares-method.html
     */
    void volatility_fitting(std::vector<double> &instvols) {

        // LMM in practice Chapter 7 page 158 step 2 - create implied caplet volatilities for a curve
        std::cout << "squared caplet implied volatilities multiplied by time to maturity" << std::endl;

        Array independentValues(capletvols.size());
        for (int i = 0; i < capletvols.size(); i++) {
            double val = tenor[i+1] * capletvols[i] * capletvols[i];
            independentValues[i] = val;
            std::cout << val << std::endl;
        }

        // initial values
        double t = 0.0; // As seen for today
        const double phi = 1.0;
        double a = 0.1; //0.112346;
        double b = 0.1; //-0.441811;
        double c = 0.1; //0.971559; //0.1;
        double d = 0.1; //1.223058; //0.1;

        // LMM in practice Chapter 7 page 158 step 3 - define fFO - Objective function with a, b, c, d parameters
        // https://www.quantlib.org/slides/dima-ql-intro-2.pdf Integral Calculation p8
        std::vector<double> squares_integrals(tenor.size(), 0.0);
        std::vector<double> difference_integrals(tenor.size(), 0.0);
        std::vector<double> fFO(independentValues.size(), 0.0);

        // instantaneous volatility curve parametric expression f(a, b, c, d , [Ti - 0])
        boost :: function < Real ( Real )> ptrF ;
        ptrF = boost :: bind (& squared_parametric_calibration , t , a ,b , c , d , _1 );

        Real absAcc = 0.00001;
        Size maxEval = 10000;
        SimpsonIntegral numInt(absAcc, maxEval);

        std::transform(tenor.begin(), tenor.end(), squares_integrals.begin(), [&](double x) {
            return numInt(ptrF, 0, x);
        });

        for (int i = 1; i < difference_integrals.size(); i++) {
            difference_integrals[i] = squares_integrals[i] - squares_integrals[i-1];
        }

        std::partial_sum(difference_integrals.begin()+1, difference_integrals.end(), fFO.begin(), std::plus<double>());

        // create corresponding curve points to be approximated
        std::vector<boost::shared_ptr<CurvePoint>> curvePoints;
        for (int i = 0; i < fFO.size(); i++) {
            //double val = parametric_calibration(0, a, b, c, d, tenor[i+1]);
            double val = fFO[i];
            curvePoints.push_back(boost::shared_ptr<CurvePoint>( new CurvePoint( val )));
        }

        // create container for function pointers for calculating rate approximations
        std::vector<boost::function<Real(const Array&)>> dependentFunctions;

        // for each curve point object, bind function pointer to operator() and add it into container
        for (unsigned int i = 0; i < curvePoints.size(); i++)
        {
            dependentFunctions.push_back(boost::bind(&CurvePoint::operator(), curvePoints[i], _1));
        }

        // perform least squares fitting and print optimized coefficients
        LeastSquares leastSquares(boost::shared_ptr<OptimizationMethod>(new LevenbergMarquardt), 10000, 1000, 1E-09, 1E-09, 1E-09);
        FittingVolConstraint parametersConstraint;
        Array initialParameters(4, 0.1);
        Array coefficients = leastSquares.Fit(independentValues, dependentFunctions, initialParameters, parametersConstraint);

#ifdef DEBUG_FITTEDVOL
        std::cout << "Results for parametric calibration coefficients" << std::endl;
        for (unsigned int i = 0; i < coefficients.size(); i++) {
            std::cout << coefficients[i] << std::endl;
        }
#endif
        a = coefficients[0];
        b = coefficients[1];
        c = coefficients[2];
        d = coefficients[3];

        instvols.resize(independentValues.size());
        for (unsigned int i = 0; i < independentValues.size(); i++) {
            instvols[i] = phi * parametric_calibration (0.0, a, b, c, d, tenor[i+1]);
        }

#ifdef DEBUG_FITTEDVOL
        // LMM in Practice Chapter 9 p159 results a,b,c,d under 9.15
        std::cout << "Results for parametric calibration" << std::endl;
        std::cout <<  "coefficients " << a  << " " << b  << " "  << c << " " << d << std::endl;
        // LMM in Practice Chapter 9 p160 table 9.16
        std::cout << "fitted volatility function" << std::endl;
        for (unsigned int i = 0; i < independentValues.size(); i++) {
            std::cout <<  "day count fraction " << tenor[i+1]  << " volatility " << instvols[i] << std::endl;
        }
#endif
    }


    /**
     * Correlation Matrix calculation
     */
    void correlation(std::vector<std::vector<double>> &rho, int size) {
        double beta = 0.1;
        rho.resize(size, std::vector<double>(size, 1.0));

        // Calculate the reminder values
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if ( (i != j) && (i > j) ){
                    double value = std::exp( -beta * (tenor[i] - tenor[j]) );
                    rho[i][j] = value;
                    rho[j][i] = value;
                }
            }
        }

#ifndef DEBUG_FITTEDVOL
        std::cout << "Correlation Matrix" << std::endl;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)
            {
                cout << rho[i][j] << " ";
            }
            std::cout << std::endl;
        }
#endif
    }

    /**
     * Run the full calibration process
     * @param instvols
     * @param rho
     * @param spot_rates
     */
    void calibrate(std::vector<double> &instvols, std::vector<std::vector<double>> &rho, std::vector<double> &spot_rates) {
        volatility_fitting(instvols);
        correlation(rho, instvols.size());
    }

private:
    double dtau;
    double expiry;
    std::vector<double> &yearFraction;
    std::vector<double> &bonds;
    std::vector<double> cplvols;
    std::vector<double> capletvols = {
        0.1641, 0.1641, 0.1641, 0.2015, 0.2189, 0.2365, 0.2550, 0.2212, 0.2255, 0.2298,
        0.2341, 0.2097, 0.2083, 0.2077, 0.2051, 0.2007, 0.1982, 0.1959, 0.1938
    };

};
