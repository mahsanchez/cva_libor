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
            volatility(instvol), rho(rho_), dtau(dtau_), expiry(expiry_)
    {
        size = instvol.size() + 1;
        initialize(spot_rates, instvol);
    }

    // LMM SDE
    void simulate(double *gaussian_rand) {
        double vol = 0.15;
        std::vector<double> dW(size, 0.0);

        // compute forward rates
        for (int n = 0; n < size-1; n++) {
            for (int i = n + 1; i < size; i++) {
                double drift = 0.0;
                for (int k = i + 1; k < size; k++) {
                    vol = volatility[k];
                    drift =+ (vol  * dtau * fwd_rates[k][n])/ (1 + dtau * fwd_rates[k][n]) * rho[i][k];
                }
                vol = volatility[n];
                double dfbar = (-drift * vol - 0.5*vol*vol) * dtau;
                dfbar += vol * gaussian_rand[n + 1] * std::sqrt(dtau);
                fwd_rates[i][n+1] = fwd_rates[i][n] * std::exp(dfbar);
            }
        }
        // compute discount factors
        for (int t = 0; t < size; t++) {
            for (int i = t + 1; i < size; i++) {
                double accuml = 1.0;
                for (int k = t; k < i; k++) {
                    accuml *= 1 / (1 + dtau * fwd_rates[k][t]);
                }
                discount_factors[i][t] = accuml;
            }
        }
#ifndef DEBUG_LMM_NUMERAIRE
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                std::cout << fwd_rates[i][j] << " ";
            }
            std::cout << std::endl;
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                std::cout << discount_factors[i][j] << " ";
            }
            std::cout << std::endl;
        }
#endif
    }

    // return forward_rates & discount_factors
    void numeraire(int index, std::vector<double> &forward_rates, std::vector<double> &discount_factor) {
        for (int t = index; t < size ; t++) {
            discount_factor[t] = discount_factors[t][index];
            forward_rates[t] = fwd_rates[t][index];
        }
    }

    inline int getSize() {
        return size;
    }

private:
    std::vector<std::vector<double>> &rho;
    std::vector<double> &volatility;
    std::vector<std::vector<double>> fwd_rates;
    std::vector<std::vector<double>> discount_factors;
    int size;
    double dtau;
    double expiry;

    // Initialize data structures
    void initialize(std::vector<double> &spot_rates, std::vector<double>& instvol) {
        // reserve memory
        fwd_rates = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
        discount_factors = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
        // Initialize fwd_rates from the spot rate curve
        for (int i = 0; i < size; i++) {
            fwd_rates[i][0] = spot_rates[i];
        }
        // initialize vol using PieceWise Constant Volatility Brigo & Mercurio p211
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
    MonteCarloSimulation(PayOff &payOff_, InterestRateModel &model_, std::vector<double> &phi_random_, int simN_) :
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
            generatePaths(&phi_random[sim * model.getSize() * model.getSize()]);

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
#ifndef DEBUG
        display_curve(exposure, "IRS Exposure");
#endif
    }

protected:
    std::vector<double>& phi_random;
    PayOff &payOff;
    InterestRateModel &model;
    std::vector<double> forward_rates;
    std::vector<double> discount_factors;
    int simN;
};



