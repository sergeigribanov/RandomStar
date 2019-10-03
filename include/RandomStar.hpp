#ifndef __RANDOMSTAR_HPP__
#define __RANDOMSTAR_HPP__
#include <vector>
#include <utility>
#include <unordered_map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <TRandom.h>
#include <CFourVector.h>
#include <MElParticle.h>

typedef std::pair<gsl_spline*, double*> InterpData;

class RandomStar {
public:
  RandomStar(const std::vector<MElParticle*>&);
  virtual ~RandomStar();
  void generate();
  const std::vector<CFourVector>& getMomenta() const;
  double getPhaseSpace() const;
  void setInitialMomentum(const CFourVector&);
private:
  double evalRoot(double, int);
  static double iBeta(double, int);
  void interpolate(int k);
  void evalPhaseSpace();
  static const std::size_t _nInterpPoints;
  TRandom _rnd;
  CFourVector _initialMomentum;
  std::vector<MElParticle*> _particles;
  int _n;
  double* _roots;
  gsl_interp_accel* _acc;
  std::unordered_map<int, InterpData> _splines;
  std::vector<CFourVector> _momenta;
  std::vector<double> _t;
  std::vector<double> _p;
  double _phspace;
  double _m_n;
};

#endif
