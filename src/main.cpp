#include <iostream>
#include <MElParticle.h>
#include "RandomStar.hpp"

int main(int argc, char* argv[]) {
  auto piPlusMeson = new MElParticle(211);
  auto piMinusMeson = new MElParticle(-211);
  auto etaMeson = new MElParticle(221);
  RandomStar star({piPlusMeson, piMinusMeson, etaMeson});
  CFourVector p(0., 0., 0., 2);
  star.setInitialMomentum(p);
  star.generate();
  std::cout << "phspace = " << star.getPhaseSpace() << std::endl;
  CFourVector psum(0., 0., 0., 0.);
  for (const auto& el : star.getMomenta()) {
    std::cout << el << std::endl;
    psum += el;
  }
  std::cout << "---- sum ---- " << std::endl;
  std::cout << psum << std::endl;
  //std::cout << "B(0.254069, 2) = " << RandomStar::iBeta(0.254069, 2) << std::endl;
  //std::cout << "x(ratio = 0.2, k = 2) = " << star.evalRoot(0.2, 2) << std::endl;
  return 0;
}
