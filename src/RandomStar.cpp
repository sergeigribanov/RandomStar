#include <iostream>
#include <algorithm>
#include <gsl/gsl_sf.h>
#include "RandomStar.hpp"

const std::size_t RandomStar::_nInterpPoints = 1000;

RandomStar::RandomStar(const std::vector<MElParticle*>& particles) :
  _particles(particles),
  _n(particles.size()),
  _roots(new double[_nInterpPoints]),
  _acc(gsl_interp_accel_alloc()),
  _momenta(_n),
  _t(_n),
  _p(_n),
  _phspace(0.) {
  double h = 1. / (_nInterpPoints - 1);
  for (std::size_t i = 0; i < _nInterpPoints; ++i) {
    _roots[i] = i * h;
  }
}

RandomStar::~RandomStar() {
  for (auto particle : _particles) {
    delete particle;
  }
  for (auto& el : _splines) {
    if (el.second.first) {
      gsl_spline_free(el.second.first);
    }
    if (el.second.second) {
      delete [] el.second.second;
    }
  }
  delete [] _roots;
  gsl_interp_accel_free(_acc);
}

void RandomStar::setInitialMomentum(const CFourVector& momentum) {
  _initialMomentum = momentum;
  _m_n = sqrt(_initialMomentum.getM2().real());
}

void RandomStar::generate() {
  double mu = 0.;
  for (const auto& el : _particles) {
    mu += el->getMass();
  }
  _t[_n - 1] = _m_n - mu;
  CFourVector p4_k = _initialMomentum;
  double m_k = _m_n;
  double m_km;
  double w_k;
  double eta_k;
  double phi_k;
  double tmp;
  int km;
  for (int k = _n; k > 1; --k) {
    km = k - 1;
    if (k > 2) {
      tmp = evalRoot(_rnd.Rndm(), k);
      _t[km - 1] = _t[km] * tmp;
    } else {
      _t[0] = 0;
    }
    tmp = _particles[km]->getMass();
    mu -= tmp;
    m_km = mu + _t[km - 1];
    w_k = 0.5 * (m_k * m_k + tmp * tmp - m_km * m_km) / m_k;
    _p[km] = sqrt(w_k * w_k - tmp * tmp);
    eta_k = 2 * _rnd.Rndm() - 1;
    phi_k = 2 * M_PI * _rnd.Rndm();
    tmp = sqrt(1 - eta_k * eta_k);
    _momenta[km][0] = _p[km] * tmp * cos(phi_k);
    _momenta[km][1] = _p[km] * tmp * sin(phi_k);
    _momenta[km][2] = _p[km] * eta_k;
    _momenta[km][3] = (p4_k.getE() * w_k +
		       p4_k.getX() * _momenta[km][0] +
		       p4_k.getY() * _momenta[km][1] +
		       p4_k.getZ() * _momenta[km][2]) / m_k;
    tmp = (_momenta[km][3].real() + w_k) / (p4_k.getE().real() + m_k);
    for (int i = 0; i < 3; ++i) {
      _momenta[km][i] += p4_k[i] * tmp;
    }
    p4_k -= _momenta[km];
    m_k = m_km;
  }
  _momenta[0] = p4_k;
  evalPhaseSpace();
}

void RandomStar::evalPhaseSpace() {
  _phspace =
    0.5 * pow(M_PI, 1.5 * (_n - 1)) *
    pow(_t[_n - 1], 1.5 * _n - 2.5) /
    _m_n / gsl_sf_gamma(1.5 * (_n - 1));
  for (int i = 1; i < _n; ++i) {
    _phspace *= _p[i] / sqrt(_t[i] - _t[i - 1]);
  }
}

double RandomStar::iBeta(double x, int k) {
  gsl_sf_result result;
  int status = gsl_sf_beta_inc_e(1.5 * (k - 1), 1.5, x, &result);
  // TODO : add exception
  return result.val;
}

double RandomStar::evalRoot(double ratio, int k) {
  if (_splines.find(k) == _splines.end()) {
    interpolate(k);
  }
  return gsl_spline_eval(_splines.at(k).first, ratio, _acc);
}

void RandomStar::interpolate(int k) {
  double fullBeta = iBeta(1, k);
  double* ratios = new double[_nInterpPoints];
  ratios[0] = 0;
  ratios[_nInterpPoints - 1] = 1.;
  for (std::size_t i = 1; i + 1 < _nInterpPoints; ++i) {
    ratios[i] = iBeta(_roots[i], k) / fullBeta;
  }
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, _nInterpPoints);
  _splines.insert(std::make_pair(k, std::make_pair(spline, ratios)));
  gsl_spline_init(spline, ratios, _roots, _nInterpPoints);
}

const std::vector<CFourVector>& RandomStar::getMomenta() const {
  return _momenta;
}

double RandomStar::getPhaseSpace() const {
  return _phspace;
}

TRandom& RandomStar::getRndGen() {
  return _rnd;
}
