#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <random>
#include "product.h"

using namespace std;

/*
 * Gaussian Random Number generators using Intel MKL
 */
class VSLRNGRandomGenerator {
public:
    VSLRNGRandomGenerator() {
        vslNewStream( &stream, VSL_BRNG_MT19937, 777 );
    }

    ~VSLRNGRandomGenerator() {
        vslDeleteStream( &stream );
    }

    void operator() (double *phi_random, int count) {
        vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, count, phi_random, 0.0, 1.0 );
    }
private:
    VSLStreamStatePtr stream;
};


/*
 * LiborMarketModelStochasticProcess
 */
class LiborMarketModel {
public:

    LiborMarketModel(std::vector<double> &spot_rates, std::vector<double>& instvol, std::vector<std::vector<double>>& rho_, double dtau_, double expiry_) :
            rho(rho_), dtau(dtau_), expiry(expiry_)
    {
        size = expiry/dtau;
         initialize(spot_rates, instvol);
    }

    // LMM SDE
    void simulate(double *gaussian_rand) {
        for (int i = 1; i < size; i++) {
            for (int t = 1; t < size; t++) {
                double drift = 0.0;
                for (int j = i+1; j < size; j++) {
                    drift += ( dtau * volatilities[j][t-1] * fwd_rates[j][t-1] )/(1 + dtau*fwd_rates[j][t-1]);
                }
                double dfbar = drift * rho[i][t-1];
                dfbar += -0.5 * volatilities[i][t-1] * volatilities[i][t-1];
                dfbar *= dtau;
                dfbar += volatilities[i][t-1] * gaussian_rand[ i * size + (t-1)] * std::sqrt(dtau);
                fwd_rates[i][t] = fwd_rates[i][t-1] * std::exp(dfbar);
            }
        }
    }

    // return forward_rates & discount_factors
    void numeraire(int index, std::vector<double> forward_rates, std::vector<double> discount_factors) {
        for (int t = index; t < size ; t++) {
            double accuml = 1.0;
            for (int k = size; k < index - 1; k--) {
                accuml *= (1 + dtau * fwd_rates[k][t]);
            }
            discount_factors[t] = 1/accuml;
            forward_rates[t] = fwd_rates[t][index];
        }
    }

    inline int getSize() {
        return size;
    }

private:
    std::vector<std::vector<double>> &rho;
    std::vector<std::vector<double>> volatilities;
    std::vector<std::vector<double>> fwd_rates;
    int size;
    double dtau;
    double expiry;

    // Initialize data structures
    void initialize(std::vector<double> &spot_rates, std::vector<double>& instvol) {
        // reserve memory
        volatilities = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
        fwd_rates = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
        // initialize vol using PieceWise Constant Volatility Brigo & Mercurio p211
        int count = size;
        for (int i = 0; i < size; i++) {
            for (int j = i; j < count; j++) {
                volatilities[j][i] = instvol[j];
            }
            count--;
        }
        // Initialize fwd_rates from the spot rate curve
        for (int i = 0; i < size; i++) {
            fwd_rates[i][0] = spot_rates[i];
        }
    }

};


/*
* Monte Carlo Simulation Engine
* Run simulation to generate the stochastics Forward Rates Risk Factors Grid using HJM Model
* Simulates forward rates using Euler-Maruyama time-stepping procedure.
* MC Apply Variance Reduction
* MC capture simulation statistics (MC Error, stddeviation, avg)
 */
template<typename InterestRateModel, typename PayOff>
class MonteCarloSimulation {
public:
    MonteCarloSimulation(PayOff &payOff_, InterestRateModel &model_, double* phi_random_, int simN_) :
            payOff(payOff_), phi_random(phi_random_), model(model_), simN(simN_)
     {
        forward_rates = std::vector<double>(model.getSize());
        discount_factors = std::vector<double>(model.getSize());
    }

    /**
     * Monte Carlo Method calculation method engine
     */
    void calculate(std::vector<std::vector<double>> &exposures, double &duration) {
        auto start = std::chrono::high_resolution_clock::now();

        // Run simulation to generate Forward Rates Risk Factors Grid using HJM Model and Musiela SDE ( CURAND )
        for (int sim = 1; sim < simN; sim++) {

            // Evolve the Forward Rates Risk Factor Simulation Path using Interest Rate Model
            generatePaths(phi_random + simN*model.getSize() * model.getSize());

            // Interest Rate Swap Mark to Market pricing the IRS across all pricing points
            pricePayOff(exposures[sim]);
        }

       // Display EE[t] curve realization for simulation sim
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        duration = elapsed.count();
    }

    /*
     * Trivial implementation to Evolve the Forward Rates Risk Factor Simulation using HJM MonteCarlo
     */
    void generatePaths(double *phi_random) {
        model.simulate(phi_random);
    }

    /*
     * Mark to Market Exposure profile calculation for Interest Rate Swap (marktomarket) pricePayOff()
     */
    void pricePayOff(std::vector<double> &exposure) {

        for (int i = 0; i < model.getSize(); i++) {
            model.numeraire(i, forward_rates, discount_factors);
            VanillaInterestRateSwapPricer pricer(payOff, forward_rates, discount_factors, model.getSize() - i);
            exposure[i] = std::max(pricer.price(), 0.0);
        }
#if DEBUG
        display_curve(exposure, "IRS Exposure");
#endif
    }

protected:
    double* phi_random;
    PayOff &payOff;
    InterestRateModel &model;
    std::vector<double> forward_rates;
    std::vector<double> discount_factors;
    int simN;
};



