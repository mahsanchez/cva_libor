#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>

#include "simulation.h"
#include "calibration_libor.h"

using namespace boost::accumulators;

// ATM Market Cap Volatility
std::vector<double> capvols_mkt = {
        0.0, 0.1641, 0.1641, 0.1641, 0.1765, 0.1889, 0.2013, 0.2137, 0.2162, 0.2186, 0.2211, 0.2235, 0.2223, 0.2212, 0.2200, 0.2188, 0.2173, 0.2158, 0.2142, 0.2127
};

// Zero Coupon Bonds
std::vector<double> zcb = {
        0.9947527, 0.9892651, 0.9834984, 0.9774658, 0.9712884, 0.9648035, 0.9580084, 0.9509789, 0.9440868, 0.9369436, 0.9295484, 0.9219838,
        0.9145031, 0.9068886, 0.8990590, 0.8911017, 0.8833709, 0.8754579, 0.8673616, 0.8581725
};


// Year Fractions
std::vector<double> yearFractions_ = {
        0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25000, 0.25278, 0.25556, 0.25556, 0.25278, 0.25278, 0.26111,
        0.25278, 0.25278, 0.25278, 0.25278, 0.25278
};

std::vector<double> yearFractions0_ = {
        0.0, 0.248804, 0.501322, 0.502435, 0.859816, 1.08345, 1.30411, 2.35441, 1.95083, 2.19285, 2.43783, 3.4725, 3.08374,  3.31253, 3.63714,
        3.97711, 4.16463, 4.41673, 4.66739, 4.81921
};


// CDS Spreads IRELAND
std::vector<double> spreads = {
        208.5, 187.5, 214, 235, 255, 272, 225, 233, 218, 215, 203, 199, 202, 196.71, 225.92, 219.69, 229.74, 232.12, 223.02, 224.45, 212,
        211.51, 206.25, 203.37, 212.94, 211.02, 210.06, 206.23, 211.49, 209.09, 204.3, 204.77, 199.98, 200.94, 202.38, 205.72, 204.76,
        210.02, 209.54, 209.54, 212.41, 213.35, 208.57, 208.56, 220.05, 226.26, 227.2, 222.89, 228.63, 231.5, 247.75, 255.37, 251.07,
        256.33, 252.01, 254.88, 246.98, 238.12, 241.95, 244.33, 252.45, 250.53, 246.71, 256.26, 255.78, 257.19, 247.63, 237.12, 234.73,
        226.36, 218, 219.9, 224.68, 221.81, 220.38, 211.77, 203.17, 206.04, 220, 218, 225, 217.5, 215, 220.5, 250.25, 260.5, 269.5, 258, 258.5,
        263, 274.5, 265.5, 268.5, 273.5, 275, 271.5, 263.75, 275, 287, 281, 271, 280.25, 284.5, 272, 275.5, 264.25, 274.25, 269.5, 264.5, 256.5,
        258, 265, 260, 262.5, 268, 272.25, 271, 274, 278.5, 278.5, 284.5, 290, 274, 264.5, 262.5, 247.75, 250.5, 248, 244.75, 246.5, 237.5, 240.5,
        236, 245.5, 237.75, 234.25, 235, 224, 215.5, 217, 217, 220.5, 208.5, 202.47, 203.43, 206.77, 205.33, 200.05, 202.91, 205.05, 222.51, 218.9,
        218.43, 221.31, 217.24, 218.67, 216.52, 216.5, 217.94, 208.37, 205.01, 200.95, 203.1, 203.81, 206.2, 204.28, 200.93, 202.36, 200.44, 197.8,
        199.23, 209.74, 211.18, 214.05, 215, 228.62, 233.63, 222.86, 218.8, 214.49, 217.36, 216.4, 213.52, 215.43, 219.49, 214.22, 218.05, 211.1,
        205.13, 207.75, 201.78, 199.39, 200.58, 199.14, 192.45, 188.38, 191.24, 192.9, 192.18, 190.26, 187.88, 186.44, 183.09, 181.18, 194.75,
        194.75, 200.75, 203, 204, 207.5, 207.25, 209, 205.75, 207.5, 203.5, 202.25, 199.5, 202, 201, 198.25, 191.5, 187.75, 188.75, 190, 193, 193.75,
        200, 200, 204, 194.5, 192.25, 189.5, 188.5, 186.5, 187.5, 193.75, 196.5, 205, 204.25, 208, 211.75, 217, 213.25, 212.5, 213.75, 211.25, 214,
        220.5, 212.5, 228, 225.75, 226.5, 233, 228.25, 225.5, 229, 229.5, 220.25, 220, 220, 223, 221, 216.5, 211.25, 199.75, 192.5, 193.94, 192.74,
        193.7, 195.6, 194.4, 195.36, 193.92, 185.79, 179.81, 162.13, 165.23, 168.81, 172.63, 168.81, 164.51, 159.01, 159.25, 154.95,
        152.08, 153.75, 153.74, 157.32, 160.66, 157.79, 152.54, 152.54, 156.84, 159.22, 162.08, 161.13, 166.38, 168.75, 73.75, 174.93, 182.8, 185.9,
        213.58, 212.15, 207.84, 214.76, 211.17, 206.62, 199.46, 196.37, 191.59, 183.95, 185.15, 172.97, 169.38, 160.08, 162.47, 159.6, 157.21, 147.91,
        145.29, 146.01, 141.95, 143.38, 136, 124.55, 119.07, 116.69, 122.17, 122.41, 122.41, 136, 136.5, 138.75, 138.5, 136.5, 135, 131.5, 130.75, 131.5,
        133.25, 133.75, 134, 136.5, 133.75, 131, 121.5, 120.5, 115.5, 114.25, 110.25, 110, 110, 110.5, 111, 111.5, 110, 112, 112, 113.25, 115.5, 120, 115.75,
        112.75, 114, 112.25, 114.5, 117, 118.25, 127, 139.5, 131, 131.75, 133, 130, 127, 132.25, 128.5, 132.25, 134.5, 136.5, 135, 138.25, 138, 134.5, 138,
        138.5, 135, 130.75, 132.5, 129.75, 125.5, 124.5, 125, 120, 119, 121.05, 122.49, 122.97, 124.88, 119.61, 119.61, 126.79, 126.78, 127.26, 127.73, 127.26,
        122.47, 118.17, 118.41, 118.65, 120.08, 123.42, 119.84, 118.64, 122.71, 124.86, 121.27, 120.79
};

struct cva_stats {
    double average;
    double median;
    double max;
    double quantile_75;
    double quantile_95;
    double quartiles;
};


/**
 * Calculate General Statistis (Mean, Percentile 95%)
 * @param stats_vector
 * @param exposures
 * @param timepoints_size
 * @param simN
 */
 /*
void report_statistics(std::vector<cva_stats>& stats_vector, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    std::vector<double> _distribution(simN);
    for (int t = 0; t < timepoints_size; t++) {
        accumulator_set<double, stats<tag::mean, tag::p_square_quantile, tag::max > > acc(quantile_probability = 0.95 );;
        for (int simulation = 0; simulation < simN; simulation++) {
            acc( exposures[simulation][t] );
        }
        stats_vector[t].average = mean(acc);
        stats_vector[t].quantile_95 = p_square_quantile(acc);
        stats_vector[t].max = boost::accumulators::max(acc);;
        acc = {};
    }
}
*/
void report_statistics(std::vector<cva_stats>& stats_vector, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    std::vector<Real> _distribution(simN);
    GeneralStatistics statistics;

    for (int t = 0; t < timepoints_size; t++) {
        for (int s = 0; s < simN; s++) {
            _distribution[s] = exposures[s][t];
        }
        statistics.addSequence(_distribution.begin(), _distribution.end());

        stats_vector[t].average = statistics.mean();
        stats_vector[t].quantile_75 = statistics.percentile(0.75);
        stats_vector[t].quantile_95 = statistics.percentile(0.95);
        stats_vector[t].max = statistics.max();

        statistics.reset();
    }
}

/*
 * Main Entry Point
 */

int main() {
    const double expiry = 5.0;
    const double dtau = 0.25;

    /*
     * Interest Rate Swap Product Definition 5Y Floating, 5Y Fixing 3M reset
     */
    std::vector<double> floating_schedule = {
            0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
    };
    std::vector<double> fixed_schedule = {
            0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
    };
    InterestRateSwap payOff(floating_schedule, floating_schedule,  fixed_schedule, 10, 0.025, expiry, dtau);

    /*
     * LMM Calibration (StrippedCapletVolatility, Fitted Instantaneous Volatility, CorrelationMatrix)
     */
    std::vector<double> vols;
    std::vector<std::vector<double>> rho;
    std::vector<double> spot_rates_;
    std::vector<double> cplvols;

    // StrippedCapletVolatility
    StrippedCapletVolatility strippedCapletVolatility(yearFractions0_,yearFractions_, zcb, capvols_mkt, expiry, dtau);
    strippedCapletVolatility.caps_strikes();
    strippedCapletVolatility.caplet_volatility_stripping(cplvols, spot_rates_);

    // Parametric Calibration to Cap Prices & Correlation Matrix
    CplVolCalibration cplVolCalibration(yearFractions, cplvols, expiry); //dtau 1Y
    cplVolCalibration.calibrate( vols, rho);

    /**
     * Monte Carlo Simulation Simulation & Exposure Profiles Generation
     */
    int simN = 3000; // Total Number of Simulations

    // Initialize the gaussians variates
    //  Increase the simulations numbers and analyze convergence and std deviation
    double duration = 0.0;
    int size = expiry/dtau;
    int count = simN * size * size;
    std::vector<double> phi_random(count, 0.0);
    VSLRNGRandomGenerator vsl_gaussian;
    vsl_gaussian(&phi_random[0], count);

    // LMM - Interest Rate  Model
    LiborMarketModel lmm(spot_rates_, vols, rho, dtau, expiry);

    // Generate one Exposure Profile by Simulation Step
    std::vector<std::vector<double>> exposures(simN, std::vector<double>(size, 0.0));

    // Monte Carlo Simulation Engine generate the Exposure IRS Grid
    MonteCarloSimulation<LiborMarketModel, InterestRateSwap> mc_engine(payOff, lmm, phi_random, simN);
    mc_engine.calculate(exposures, duration);

#ifdef DEBUG_EXPOSURE_CVA
    std::cout << "Exposures Profile" << std::endl;
    for (int i = 0; i < exposures.size(); i++) {
        for (int j = 0; j < exposures[j].size(); j++) {
            std::cout << exposures[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif

    /*
     * Counter Party Credit Risk Analytics & Expected Exposure Profile
     */
    double recovery = 0.04;      // Recovery Rates
    size = tenor.size() - 1;

    // Expected Exposure Profile
    std::vector<cva_stats> stats_vector;
    stats_vector.resize(size);

    // Calculate Expected Exposure (EE)  EPE(t) = 𝔼 [max(V , 0)|Ft]
    // Use reduction to generate the expectation on the distribution across timepoints the expected exposure profile for the IRS
    // Calculate Statistics max, median, quartiles, 97.5th percentile on exposures
    // Calculate Potential Future Exposure (PFE) at the 97.5th percentile and media of the Positive EE
    report_statistics(stats_vector, exposures, size, simN);

    std::vector<double> expected_exposure(tenor.size(), 0.0);
    std::transform(stats_vector.begin(), stats_vector.end(), expected_exposure.begin(), [](cva_stats s) {
        return s.average;
    });

    // Report Expected Exposure Profile Curve
#ifdef DEBUG_EXPOSURE_CVA
    std::cout << "Tenors" << std::endl;
    for (int j = 0; j < tenor.size(); j++) {
        std::cout << tenor[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "Expected Exposure EE[t]" << std::endl;
    for (int j = 0; j < expected_exposure.size(); j++) {
        std::cout << expected_exposure[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "Potential Future Exposure (0.95) PFE[t]" << std::endl;
    for (int j = 0; j < stats_vector.size(); j++) {
        std::cout << stats_vector[j].quantile_95 << " ";
    }
    std::cout << std::endl;
#endif

    /*
     * Counter Party Credit Risk Analytics & Credit Value Adjustment
     */

    InterpolatedFICurve zcbCurve(tenor, zcb);
    InterpolatedFICurve spotRatesCurve(tenor, spot_rates_);
    InterpolatedSpreadsCurve spreadsCurve(tenor, spreads, 0.0001);

    // Survival Probability Bootstrapping
    SurvivalProbCurve survivalProbCurve(tenor, spreadsCurve, zcbCurve, recovery);

#ifdef DEBUG_CDS_CVA
    std::cout << "CDS - Survival Probabilities" << std::endl;
    for (int i = 0; i < tenor.size(); i++) {
        std::cout << tenor[i] << " " << survivalProbCurve(tenor[i]) << std::endl;
    }
    std::cout << std::endl;
#endif

    // Calculate The Unilateral CVA - Credit Value Adjustment Metrics Calculation.
    // For two conterparties A - B. Counterparty A want to know how much is loss in a contract due to B defaulting

    /*
    expected_exposure = {0.0, 0.253899346, 0.262649146,0.264752078 , 0.269196504, 0.26545791, 0.258220705,	0.251405001, 0.244940384, 0.232968488, 0.21709824,
                         0.201386603, 0.184246938, 0.164233229, 0.142736561, 0.121057271, 0.098013021, 0.07231877, 0.042947986, 0.0 };
    */

    InterpolatedFICurve expectedExposureCurve(tenor, expected_exposure);

    /*
     * CVA Calculation
     * CVA =  E [ (1 - R) [ DF[t] * EE[t] * dPD[t] ] ]
     */
    std::vector<double> defaultProbabilities;
    std::transform(tenor.begin(), tenor.end(), std::back_inserter(defaultProbabilities), [&](double t) {
        return 1.0 - survivalProbCurve(t);
    });

    double cva = 0.0;
    for (int i = 1; i < tenor.size(); i++) {
        cva += (1 - recovery) * zcbCurve(tenor[i]) * expectedExposureCurve(tenor[i]) * (defaultProbabilities[i] - defaultProbabilities[i-1]);
    }

#ifdef DEBUG_EXPOSURE_CVA
    std::cout << "Credit value Adjustment " << std::endl;
    std::cout << std::setprecision(6)<< std::fixed << cva << " " << simN << " " << duration << std::endl;
#endif

    exit(0);
}