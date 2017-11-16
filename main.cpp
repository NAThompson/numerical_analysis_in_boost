#include <iostream>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using boost::math::cubic_b_spline;
using boost::math::barycentric_rational;
using boost::math::chebyshev_transform;
using boost::math::quadrature::trapezoidal;
using boost::math::quadrature::tanh_sinh;
using boost::math::quadrature::gauss_kronrod;
using boost::math::constants::two_pi;
using boost::math::constants::half_pi;

int main()
{
    std::vector<double> v(500);
    for (size_t i = 0; i < v.size(); ++i)
    {
        v[i] = std::sin(i*M_PI/v.size());
    }
    cubic_b_spline<double> s(v.begin(), v.end(), 0, 1.0);

    std::cout << "s(" << 0 << ") = " << s(0) << std::endl;
    std::cout << "s'(" << 0 << ") = " << s.prime(0) << std::endl;

    std::vector<double> x(500);
    std::vector<double> y(500);
    for (size_t i = 0; i < v.size(); ++i)
    {
        x[i] = i*M_PI/v.size();
        y[i] = std::sin(i*M_PI/v.size());
    }

    barycentric_rational<double> r(x.data(), y.data(), x.size(), 4);

    std::cout << "r(" << 0.5 << ") = " << r(0.5) << std::endl;


    auto f = [](double x) { return sin(x); };

    chebyshev_transform<double> cheb(f, 0.0, M_PI);
    std::cout << "cheb(" << 0.7 << ") = " << cheb(0.7) << std::endl;
    std::cout << "cheb'(" << 0.7 << ") = " << cheb.prime(0.7) << std::endl;
    std::cout << "I[cheb] = " << cheb.integrate() << std::endl;


    auto g = [](double x) { return 1/(5 - 4*cos(x)); };
    double I = trapezoidal(g, 0.0, two_pi<double>());
    std::cout << "I = " << I << std::endl;

    tanh_sinh<double> integrator;
    auto f2 = [&](double x)->double { return 1/(1 + pow(tan(x), 3)); };
    double Q = integrator.integrate(f2, (double) 0, half_pi<double>());
    std::cout << "Q = " << Q << std::endl;
    gauss_kronrod<double, 15> gk15;
    Q = gk15.integrate(f2, 0.0, half_pi<double>());
    std::cout << "Q = " << Q << " by Gauss-Kronrod\n";
}
