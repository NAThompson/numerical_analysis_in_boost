#include <iostream>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

using boost::math::cubic_b_spline;
using boost::math::barycentric_rational;
using boost::math::chebyshev_transform;
using boost::math::quadrature::trapezoidal;

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
    double I = trapezoidal(g, 0.0, boost::math::constants::two_pi<double>());
    std::cout << "I = " << I << std::endl;

}
