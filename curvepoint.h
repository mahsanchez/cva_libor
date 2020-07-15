#include <ql/quantlib.hpp>
using namespace QuantLib;
//
// class modeling rate approximation for a single
// curve point using a set of given parameters
class CurvePoint
{
public:
    CurvePoint(Real t) : t(t) { }

    Real operator()(const Array& coefficients){
        Real approximation = 0.0;
        for (unsigned int i = 0; i < coefficients.size(); i++)
        {
            approximation += coefficients[i] * pow(t, i);
        }
        return approximation;
    }
private:
    Real t;
};