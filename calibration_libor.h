#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>

#include "mkl_lapacke.h"

#include <ql/quantlib.hpp>
#include <ql/math/solvers1d/bisection.hpp>

using namespace QuantLib;

// number of days since (Ti - T0)/360
std::vector<double> tenor = {
    //3M, 6M, 9M, 1Y, 1Y3M, 1Y6M, 1Y9M, 2Y, 2Y3M, 2Y6M, 2Y9M, 3Y, 3Y3M, 3Y6M, 3Y9M, 4Y, 4Y3M, 4Y6M, 4Y9M, 5Y
    0.2514, 0.5028, 0.7583, 1.0139, 1.2639, 1.5167, 1.7722, 2.0278, 2.2778, 2.5306, 2.7861, 3.0417, 3.2944, 3.5472, 3.8083, 4.0611, 4.3139, 4.5667, 4.8194, 5.0722
};

std::vector<double> yearFractions = {
     0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25278, 0.25278, 0.26111,
     0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278, 0.25278
};

/*
 * Price a CAPLET using Black Formula
 */
double BlackIRCaplet(double zcb, double yearFraction0, double yearFraction, double strike, double x, double fwd_rate, double vol) {
    double d1 = std::log(fwd_rate/x) + 0.5*vol*vol * yearFraction0;
    d1 /= vol * std::sqrt(yearFraction0);
    double d2 = d1 - vol * std::sqrt(yearFraction0);
    double Nd1 = 0.0;
    double Nd2 = 0.0;
    vdCdfNorm( 1, &d1, &Nd1 );
    vdCdfNorm( 1, &d2, &Nd2 );
    double result = fwd_rate * Nd1 - strike * Nd2;
    result *= zcb;
    result *= yearFraction;
    return result;
}

/**
 * Implied Caplet Volatility Problem LMM in Practice Algorithm 7.1
 */
double impliedVolProblem(double zcb, double yearFraction0, double yearFraction, double strike, double x, double fwd_rate, double vol, double cap, double cplet_price) {
    return (cap - cplet_price - BlackIRCaplet(zcb, yearFraction0, yearFraction, strike, x, fwd_rate, vol));
}

/**
 * Stripped Caplet Volatility
 */

class StrippedCapletVolatility {
public:
    StrippedCapletVolatility(std::vector<double> &yearFraction0_, std::vector<double> &yearFraction_, std::vector<double> &zcb_, std::vector<double> &capvols_, double expiry_, double dtau_) :
        yearFraction0(yearFraction0_), yearFraction(yearFraction_), zcb(zcb_), capvols_mkt(capvols_), expiry(expiry_), dtau(dtau_) {
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
    void caplet_volatility_stripping(std::vector<double> &cplvols, std::vector<double> &spot_rates) {
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
            acc[i] = BlackIRCaplet(zcb[i], yearFraction0[i], yearFraction[i], strikes[i], strikes[i], spot_rates[i], capvols_mkt[i]);
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

        // LMM in Practice chapter 7 p76 Algorithm 7.1
        cplvols.resize(size);
        cplvols[1] = 0.1641;     // 6M

        for (int i = 2; i < cplvols.size(); i++) {
            double cplet_price = 0.0; // caplet_price
            for (int j = 1; j < i; j++) {
                cplet_price += BlackIRCaplet(zcb[j], yearFraction0[j], yearFraction[j], strikes[j], strikes[j], spot_rates[j], cplvols[j]);
            }
#ifdef DEBUG_IMPLIEDCAPLETVOL
            std::cout << "caps " << caps_quotes[i] << " caplet price " << cplet_price << std::endl;
#endif
            //calculate implied caplet volatility
            boost :: function < Real ( Volatility )> impliedVolFunc ;   // setup a boost function
            // bind the boost function to all market parameters , keep vol as variant
            impliedVolFunc = boost :: bind (& impliedVolProblem , zcb[i], yearFraction0[i], yearFraction[i], strikes[i], strikes[i], spot_rates[i], _1, caps_quotes[i], cplet_price);
            // root finding using interval bisection algorithm
            cplvols[i] = root_solver.solve(impliedVolFunc, 1E-6, 0.12, 0.10, 0.35);
        }

#ifdef DEBUG_IMPLIEDCAPLETVOL
        // LMM in Practice chapter 7 p77 Table 7.5
        std::cout << "Caplet Volatilites Stripped from cap volatilities (cap vol, refernce implied capvol, implied cpvol)" << std::endl;
        for (int i = 1; i < cplvols.size(); i++) {
            std::cout << tenor[i]  << " " << capvols_mkt[i]  << " " << capletvols[i] << " " << cplvols[i] << std::endl;
        }
#endif
    }

private:
    double dtau;
    double expiry;
    std::vector<double> &yearFraction0;
    std::vector<double> &yearFraction;
    std::vector<double> &zcb;
    std::vector<double> &capvols_mkt;
    std::vector<double> strikes;
    std::vector<double> capletvols = {
            0.0, 0.1641, 0.1641, 0.1641, 0.2015, 0.2189, 0.2365, 0.2550, 0.2212, 0.2255, 0.2298,
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
 * Least Squares Function
 */
class LeastSquares : public CostFunction {
public:
    LeastSquares(Array& independentValues_, std::vector<double> &tenor_) : independentValues(independentValues_), tenor(tenor_){
    }

    Real value(const Array& x) const {
        Array differences = values(x);
        Real sumOfSquaredDifferences = 0.0;
        for (unsigned int i = 0; i < differences.size(); i++)
        {
            sumOfSquaredDifferences += differences[i] * differences[i];
        }

        std::cout << "Cost Function " << std::sqrt(sumOfSquaredDifferences) << " " << x << std::endl;

        return std::sqrt(sumOfSquaredDifferences);
    }

    // LMM in practice Chapter 7 page 158 step 3 - define fFO - Objective function with a, b, c, d parameters
    // https://www.quantlib.org/slides/dima-ql-intro-2.pdf Integral Calculation p8
    // calculate differences between all observed and estimated values using
    // function pointers to calculate estimated value using a set of given parameters
    Disposable<Array> values(const Array& x) const
    {
        // fFO[i] Vector calculation
        // instantaneous volatility curve parametric expression f(t, a, b, c, d , [Ti - t])
        std::vector<double> squares_integrals(tenor.size(), 0.0);
        std::vector<double> difference_integrals(tenor.size(), 0.0);
        std::vector<double> fFO(independentValues.size(), 0.0);

        boost :: function < Real ( Real )> ptrF ;
        ptrF = boost :: bind (& squared_parametric_calibration , 0.0 , x[0], x[1], x[2], x[3], _1 );

        Real absAcc = 0.00001;
        Size maxEval = 10000;
        SimpsonIntegral numInt(absAcc, maxEval);

        std::transform(tenor.begin(), tenor.end(), squares_integrals.begin(), [&](double value) {
            return numInt(ptrF, 0, value);
        });

        for (int i = 1; i < difference_integrals.size(); i++) {
            difference_integrals[i] = squares_integrals[i] - squares_integrals[i-1];
        }

        std::partial_sum(difference_integrals.begin()+1, difference_integrals.end(), fFO.begin());

        // Least Squares Differences
        Array differences(independentValues.size());
        for (unsigned int i = 0; i < independentValues.size(); i++)
        {
            differences[i] = independentValues[i] - fFO[i];
        }

        return differences;
    }


    virtual Real valueAndGradient(Array &grad, const Array& x) const {
        gradient(grad, x);
        return value(x);
    }

private:
    Array& independentValues;
    std::vector<double> &tenor;
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
            if ( (params[0] + params[1]) <= 0.0) return false;
            if (params[3] <= 0.0) return false;
            return true;
        }
        Array upperBound(const Array& params) const {
            return Array(params.size(), std::numeric_limits < Array::value_type > ::max()); //
        }
        Array lowerBound(const Array& params) const {
            Array lower_bound(params.size(), -2.0);
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
    CplVolCalibration(std::vector<double> &yearFraction_, std::vector<double> &capletvols_, double expiry_) :
       yearFraction(yearFraction_), cplvols(capletvols_), expiry(expiry_), dtau(0.25) {}

    /*
     * Volatility Fitting - solve optimization problem  argmin sum[cplvols - fitvols]^2 using least square methods
     * LMM in practice Chapter 9
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

        //instantiate Constraint
        FittingVolConstraint parametricConstraint;

        Size maxIterations = 100000; // end search if after 10000 iterations if no solution
        Size minStatIterations = 10; // don't spend more than 10 iterations in a single point
        Real rootEpsilon = 1e-5; // end search if absolute difference of current and last fraction value is under a threshold
        Real functionEpsilon = 1e-5; // end search if absolute difference of current and last function value is under a threshold
        Real gradientEpsilon = 1e-5; //0.00001; // end search if absolute difference of norm of current and last gradient is bellow epsilum

        EndCriteria endCriteria(maxIterations, minStatIterations, rootEpsilon, functionEpsilon, gradientEpsilon);
        Array InitVal(4, 0.1);

        LeastSquares leastSquares(independentValues, tenor);
        Problem getValues(leastSquares, parametricConstraint, InitVal);

        // Methods
        SteepestDescent sd;
        ConjugateGradient cg;
        BFGS sol;
        LevenbergMarquardt lm;

        // if the algorithm is able to locate an optimal solution it will stop searching at the stationary point in the search space
        EndCriteria::Type solution = sd.minimize(getValues, endCriteria);

        Array coefficients = getValues.currentValue();

        // LMM in Practice Chapter 9 p159 results a,b,c,d under 9.15
#ifdef DEBUG_FITTEDVOL
        std::cout << "Results for parametric calibration coefficients" << std::endl;
        std::cout << "Solution type: %" << std::endl;
        std::cout << "Root: %" << coefficients << std::endl;
        std::cout << "Reference Min F Value: 0.0436646" << " Actual Min F Value: %" << getValues.functionValue() << std::endl;
#endif

        double a = coefficients[0];  //0.112346;
        double b = coefficients[1];  //-0.441811;
        double c = coefficients[2];  //0.971559;
        double d = coefficients[3];  //1.223058;

        std::vector<double> reference_instvols(independentValues.size());
        instvols.resize(independentValues.size());
        for (unsigned int i = 0; i < independentValues.size(); i++) {
            instvols[i] = phi * parametric_calibration (0.0, a, b, c, d, tenor[i+1]);
            reference_instvols[i] = phi * parametric_calibration (0.0, 0.112346, -0.441811, 0.971559, 1.223058, tenor[i+1]);
        }

#ifdef DEBUG_FITTEDVOL
        std::cout <<  "reference coefficients " << 0.112346  << " " << -0.441811  << " "  << 0.971559 << " " << 1.223058 << std::endl;
        std::cout <<  "actual coefficients " << a  << " " << b  << " "  << c << " " << d << std::endl;
        // LMM in Practice Chapter 9 p160 table 9.16
        std::cout << "fitted volatility function" << std::endl;
        for (unsigned int i = 0; i < instvols.size(); i++) {
            std::cout <<  "tenor " << tenor[i+1]  << " reference instantaneous volatility " << reference_instvols[i] <<  " actual instantaneous volatility " << instvols[i] << std::endl;
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

#ifdef DEBUG_FITTEDVOL
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
     */
    void calibrate(std::vector<double> &instvols, std::vector<std::vector<double>> &rho) {
        volatility_fitting(instvols);
        correlation(rho, instvols.size());
    }

private:
    double dtau;
    double expiry;
    std::vector<double> &yearFraction;
    std::vector<double> cplvols;
    std::vector<double> capletvols = {
        0.1641, 0.1641, 0.1641, 0.2015, 0.2189, 0.2365, 0.2550, 0.2212, 0.2255, 0.2298,
        0.2341, 0.2097, 0.2083, 0.2077, 0.2051, 0.2007, 0.1982, 0.1959, 0.1938
    };

};
