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

const int tenor_size = 51;

// First row is the last observed forward curve (BOE data)  //  3-dimensional Normal random vector in columns BC, BD, BE (on far right)
std::vector<double> spot_rates = {
        0.046138361,0.045251174,0.042915805,0.04283311,0.043497719,0.044053792,0.044439518,0.044708496,0.04490347,0.045056615,0.045184474,0.045294052,0.045386152,0.045458337,0.045507803,0.045534188,0.045541867,0.045534237,0.045513128,0.045477583,0.04542292,0.045344477,0.04523777,0.045097856,0.044925591,0.04472353,0.044494505,0.044242804,0.043973184,0.043690404,0.043399223,0.043104398,0.042810688,0.042522852,0.042244909,0.041978295,0.041723875,0.041482518,0.04125509,0.041042459,0.040845492,0.040665047,0.040501255,0.040353009,0.040219084,0.040098253,0.039989288,0.039890964,0.039802053,0.039721437,0.03964844
};


// CDS Spreads IRELAND
std::vector<double> spreads = {
        208.5, 187.5, 214, 235, 255, 272, 225, 233, 218, 215, 203, 199, 202, 196.71, 225.92, 219.69, 229.74, 232.12, 223.02, 224.45, 212, 211.51, 206.25, 203.37, 212.94, 211.02, 210.06, 206.23, 211.49, 209.09, 204.3, 204.77,
        199.98, 200.94, 202.38, 205.72, 204.76, 210.02, 209.54, 209.54, 212.41, 213.35, 208.57, 208.56, 220.05, 226.26, 227.2, 222.89, 228.63, 231.5, 247.75, 255.37, 251.07,
        256.33, 252.01, 254.88, 246.98, 238.12, 241.95, 244.33, 252.45, 250.53, 246.71, 256.26, 255.78, 257.19, 247.63, 237.12, 234.73, 226.36, 218, 219.9, 224.68,
        221.81, 220.38, 211.77, 203.17, 206.04, 220, 218, 225, 217.5, 215, 220.5, 250.25, 260.5, 269.5, 258, 258.5, 263, 274.5, 265.5, 268.5, 273.5, 275, 271.5, 263.75, 275, 287, 281, 271, 280.25, 284.5, 272, 275.5, 264.25, 274.25,
        269.5, 264.5, 256.5, 258, 265, 260, 262.5, 268, 272.25, 271, 274, 278.5, 278.5, 284.5, 290, 274, 264.5, 262.5, 247.75, 250.5, 248, 244.75, 246.5, 237.5, 240.5, 236,
        245.5, 237.75, 234.25, 235, 224, 215.5, 217, 217, 220.5, 208.5, 202.47, 203.43, 206.77, 205.33, 200.05, 202.91, 205.05, 222.51, 218.9, 218.43, 221.31, 217.24, 218.67,
        216.52, 216.5, 217.94, 208.37, 205.01, 200.95, 203.1, 203.81, 206.2, 204.28, 200.93, 202.36, 200.44, 197.8, 199.23, 209.74, 211.18, 214.05, 215, 228.62, 233.63, 222.86, 218.8, 214.49, 217.36, 216.4, 213.52, 215.43,
        219.49, 214.22, 218.05, 211.1, 205.13, 207.75, 201.78, 199.39, 200.58, 199.14, 192.45, 188.38, 191.24, 192.9, 192.18, 190.26, 187.88, 186.44, 183.09, 181.18, 194.75,
        194.75, 200.75, 203, 204, 207.5, 207.25, 209, 205.75, 207.5, 203.5, 202.25, 199.5, 202, 201, 198.25, 191.5, 187.75, 188.75, 190, 193, 193.75, 200, 200, 204, 194.5,
        192.25, 189.5, 188.5, 186.5, 187.5, 193.75, 196.5, 205, 204.25, 208, 211.75, 217, 213.25, 212.5, 213.75, 211.25, 214, 220.5, 212.5, 228, 225.75, 226.5, 233, 228.25, 225.5, 229, 229.5, 220.25, 220, 220, 223, 221, 216.5,
        211.25, 199.75, 192.5, 193.94, 192.74, 193.7, 195.6, 194.4, 195.36, 193.92, 185.79, 179.81, 162.13, 165.23, 168.81, 172.63, 168.81, 164.51, 159.01, 159.25, 154.95,
        152.08, 153.75, 153.74, 157.32, 160.66, 157.79, 152.54, 152.54, 156.84, 159.22, 162.08, 161.13, 166.38, 168.75, 73.75, 174.93, 182.8, 185.9, 213.58, 212.15,
        207.84, 214.76, 211.17, 206.62, 199.46, 196.37, 191.59, 183.95, 185.15, 172.97, 169.38, 160.08, 162.47, 159.6, 157.21, 147.91, 145.29, 146.01, 141.95, 143.38, 136, 124.55, 119.07, 116.69, 122.17, 122.41, 122.41, 136, 136.5,
        138.75, 138.5, 136.5, 135, 131.5, 130.75, 131.5, 133.25, 133.75, 134, 136.5, 133.75, 131, 121.5, 120.5, 115.5, 114.25, 110.25, 110, 110, 110.5, 111, 111.5, 110,
        112, 112, 113.25, 115.5, 120, 115.75, 112.75, 114, 112.25, 114.5, 117, 118.25, 127, 139.5, 131, 131.75, 133, 130, 127, 132.25, 128.5, 132.25, 134.5, 136.5, 135, 138.25,
        138, 134.5, 138, 138.5, 135, 130.75, 132.5, 129.75, 125.5, 124.5, 125, 120, 119, 121.05, 122.49, 122.97, 124.88, 119.61, 119.61, 126.79, 126.78, 127.26, 127.73, 127.26, 122.47, 118.17, 118.41, 118.65, 120.08, 123.42, 119.84,
        118.64, 122.71, 124.86, 121.27, 120.79
};

struct cva_stats {
    double average;
    double median;
    double max;
    double quantile_75;
    double quartiles;
};


void report_statistics(std::vector<cva_stats>& stats_vector, std::vector<std::vector<double>>& exposures, int timepoints_size, int simN) {
    //std::vector<double> mc_datapoints(simN);
    for (int t = 0; t < timepoints_size; t++) {
        accumulator_set<double, stats<tag::mean, tag::p_square_quantile, tag::max > > acc(quantile_probability = 0.75 );;
        for (int simulation = 0; simulation < simN; simulation++) {
            //mc_datapoints[simulation] = exposures[t][simulation];
            acc( exposures[t][simulation] );
        }
        stats_vector[t].average = mean(acc);
        stats_vector[t].quantile_75 = p_square_quantile(acc);
        stats_vector[t].max = boost::accumulators::max(acc);;
        acc = {};
#ifndef DEBUG_STATISTICS
#endif
    }
    int a = 0;
}

/**
 *
 */

int main() {
    const double expiry = 5.0;
    const double dtau = 0.25;

    // Interest Rate Swap Product Definition 5Y Floating, 5Y Fixing 3M reset
    std::vector<double> floating_schedule = {
            0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
    };
    std::vector<double> fixed_schedule = {
            0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0
    };
    InterestRateSwap payOff(floating_schedule, floating_schedule,  fixed_schedule, 10, 0.025, expiry, dtau);

    //  Increase the simulations numbers and analyze convergence and std deviation
    int simN = 3000;
    double duration = 0.0;

    /*
     * LMM Calibration (VolatilityCurve, CorrelationMatrix)
     */
    std::vector<double> vols;
    std::vector<std::vector<double>> rho;
    std::vector<double> spot_rates_;

    // Instantaneous Volatility Calibration and Correlation Matrix
    CplVolCalibration cplVolCalibration(yearFractions, bonds, capvols, expiry); //dtau 1Y
    cplVolCalibration.calibrate( vols, rho, spot_rates_);

    // Initialize the gaussians variates
    int size = expiry/dtau;
    int count = simN * size * size;
    std::vector<double> phi_random(count, 0.0);
    VSLRNGRandomGenerator vsl_gaussian;
    vsl_gaussian(&phi_random[0], count);

    // LMM - Interest Rate  Model
    LiborMarketModel lmm(spot_rates, vols, rho, dtau, expiry);

    // Generate one Exposure Profile by Simulation Step
    std::vector<std::vector<double>> exposures(simN, std::vector<double>(size, 0.0));

    // Monte Carlo Simulation Engine generate the Exposure IRS Grid
    MonteCarloSimulation<LiborMarketModel, InterestRateSwap> mc_engine(payOff, lmm, phi_random, simN);
    mc_engine.calculate(exposures, duration);

#ifdef DEBUG_EXPOSURE_CVA
    for (int i = 0; i < exposures.size(); i++) {
        for (int j = 0; j < exposures[j].size(); j++) {
            std::cout << exposures[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif


    /*
     * Counter Party Credit Risk Analytics
     */

    // Recovery Rates
    double recovery = 0.04;

    // Expected Exposure Profile
    std::vector<cva_stats> stats_vector;
    stats_vector.resize(size);

    // Calculate Expected Exposure (EE)  EPE(t) = 𝔼 [max(V , 0)|Ft]
    // Use reduction to generate the expectation on the distribution across timepoints the expected exposure profile for the IRS
    // Calculate Statistics max, median, quartiles, 97.5th percentile on exposures
    report_statistics(stats_vector, exposures, size, simN);

    std::vector<double> expected_exposure(size, 0.0);
    std::transform(stats_vector.begin(), stats_vector.begin() + size, expected_exposure.begin(), [](cva_stats s) {
        return s.average;
    });

    // Report Expected Exposure Profile Curve
    display_curve(expected_exposure, "Expected Exposure");

    // Calculate Potential Future Exposure (PFE) at the 97.5th percentile and media of the Positive EE

    // Discounts  or Zero Coupon Bond bootstrapped from the Spot Rates
    SpotRateYieldCurveTermStructure yieldCurve(spot_rates, expiry, dtau);

    // Survival Probability Bootstrapping
    SurvivalProbabilityTermStructure survivalProbabilityCurve(tenors, spreads, yieldCurve, recovery, tenor_size); // dtau

    // Calculate The Unilateral CVA - Credit Value Adjustment Metrics Calculation.
    // For two conterparties A - B. Counterparty A want to know how much is loss in a contract due to B defaulting
    ExpectedExposureTermStructure expectedExposureCurve(floating_schedule,expected_exposure, expiry);
    double cva = calculate_cva(recovery, yieldCurve, expected_exposure, survivalProbabilityCurve, floating_schedule, expiry);

    //std::cout << std::setprecision(6)<< std::fixed << cva << " " << simN << " " << duration << std::endl;


    exit(0);
}