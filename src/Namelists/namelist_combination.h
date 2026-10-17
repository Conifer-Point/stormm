// -*-c++-*-
#ifndef STORMM_NAMELIST_COMBINATION_H
#define STORMM_NAMELIST_COMBINATION_H

#include "copyright.h"
#include "Potential/pme_util.h"
#include "nml_dynamics.h"
#include "nml_minimize.h"
#include "nml_pppm.h"

namespace stormm {
namespace namelist {

/// \brief Determine the value of a cutoff based on information given in either of two namelists.
///        STORMM will allow users to specify the cutoffs for electrostatic or Lennard-Jones
///        interactions in &dynamics, &minimize, or &pppm namelists, whichever is most intuitive.
///        However, there is a priority system in the event that specifications conflict.  A
///        defined value will always superceded a default value.  Information in the &pppm namelist
///        control block, assumed to be the domain of more advanced users, will supercede any
///        equivalent information given in the &dynamics or &minimize namelists.
///
/// Overloaded:
///   - Accept one &pppm namelist.  The first namelist must describe electrostatics,
///   - Accept two &pppm namelists.  One of the namelists must describe electrostatics while the
///     other describes van-der Waals interactions.  
///
/// \param foocon   Parameters from a &dynamics (molecular dynamics) or &minimize (energy
///                 minimization) control block
/// \param pppmcon  Parameters from a &pppm (Particle-Particle, Particle-Mesh) control block
/// \param theme    Indicate whether to arbitrate the non-bonded cutoff for electrostatic or
///                 van-der Waals interactions.  A value of ALL is an error in this context.
/// \{
template <typename T>
double arbitrateCutoff(const T &foocon, const PPPMControls &pmecon, NonbondedTheme theme);

template <typename T>
double arbitrateCutoff(const T &foocon, const PPPMControls &pmecon_a, const PPPMControls &pmecon_b,
                       NonbondedTheme theme);

template <typename T>
double arbitrateCutoff(const T &foocon, const std::vector<PPPMControls> &pmeconv,
                       NonbondedTheme theme);
/// \}

/// \brief Determine the Ewald coefficient for particle-mesh electrostatics based on information
///        given in either of two namelists.  The particle-particle interaction table derives its
///        Ewald coefficient from the electrostatic cutoff (see arbitrateCutoff(), above) and the
///        direct sum tolerance, and the reciprocal-space convolution must use the same value.
///        The Ewald coefficient stored in the &pppm control block derives from the cutoff given
///        in that block, or its default, and differs whenever the arbitrated cutoff came from the
///        &dynamics or &minimize namelist and is not the same.
///
/// \param foocon  Parameters from a &dynamics or &minimize control block
/// \param pmecon  Parameters from a &pppm control block, including the direct sum tolerance
template <typename T>
double arbitrateEwaldCoefficient(const T &foocon, const PPPMControls &pmecon);

} // namespace namelist
} // namespace stormm

#include "namelist_combination.tpp"

#endif
