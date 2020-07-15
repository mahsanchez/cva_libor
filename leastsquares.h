#include <ql/quantlib.hpp>
using namespace QuantLib;
//
class LeastSquares
{
public:
 LeastSquares(boost::shared_ptr<OptimizationMethod> solver, Size maxIterations, Size maxStationaryStateIterations, 
  Real rootEpsilon, Real functionEpsilon, Real gradientNormEpsilon) : solver(solver), maxIterations(maxIterations), maxStationaryStateIterations(maxStationaryStateIterations),
           rootEpsilon(rootEpsilon), functionEpsilon(functionEpsilon), gradientNormEpsilon(gradientNormEpsilon) { }
 //
 Array Fit(const Array& independentValues, std::vector<boost::function<Real(const Array&)>> dependentFunctions, 
  const Array& initialParameters, Constraint parametersConstraint) {
     ObjectiveFunction objectiveFunction(dependentFunctions, independentValues);
     EndCriteria endCriteria(maxIterations, maxStationaryStateIterations, rootEpsilon, functionEpsilon, gradientNormEpsilon);
     Problem optimizationProblem(objectiveFunction, parametersConstraint, initialParameters);
     //LevenbergMarquardt solver;
     EndCriteria::Type solution = solver->minimize(optimizationProblem, endCriteria);
     return optimizationProblem.currentValue();
 }

private:
 boost::shared_ptr<OptimizationMethod> solver;
 Size maxIterations; 
 Size maxStationaryStateIterations;
 Real rootEpsilon; 
 Real functionEpsilon;
 Real gradientNormEpsilon;
 //
 // cost function implementation as private nested class 
     class ObjectiveFunction : public CostFunction
     {
     private:
      std::vector<boost::function<Real(const Array&)>> dependentFunctions;
      Array independentValues;
     public:
      ObjectiveFunction(std::vector<boost::function<Real(const Array&)>> dependentFunctions, const Array& independentValues) : dependentFunctions(dependentFunctions), independentValues(independentValues) { }

      // calculate squared sum of differences between observed and estimated values
      Real value(const Array& x) const {
          Array differences = values(x);
          Real sumOfSquaredDifferences = 0.0;
          for (unsigned int i = 0; i < differences.size(); i++)
          {
              sumOfSquaredDifferences += differences[i] * differences[i];
          }
          return sumOfSquaredDifferences;
      }

      Disposable<Array> values(const Array& x) const
      {
          // calculate differences between all observed and estimated values using
          // function pointers to calculate estimated value using a set of given parameters
          Array differences(dependentFunctions.size());
          for (unsigned int i = 0; i < dependentFunctions.size(); i++)
          {
              differences[i] = dependentFunctions[i](x) - independentValues[i];
          }
          return differences;
      }
     };

};