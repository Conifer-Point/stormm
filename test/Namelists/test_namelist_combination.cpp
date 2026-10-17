#include <cmath>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Namelists/namelist_combination.h"
#include "../../src/Namelists/namelist_emulator.h"
#include "../../src/Namelists/namelist_enumerators.h"
#include "../../src/Namelists/nml_dynamics.h"
#include "../../src/Namelists/nml_pppm.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::ExceptionResponse;
using stormm::constants::tiny;
using stormm::energy::default_pme_cutoff;
using stormm::energy::ewaldCoefficient;
using stormm::energy::NonbondedTheme;
using stormm::parse::NumberFormat;
using stormm::parse::realToString;
using stormm::parse::TextFile;
using stormm::parse::TextOrigin;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::namelist;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Ewald coefficient arbitration between &dynamics and &pppm");

  // A cutoff given in the &dynamics namelist takes precedence over the &pppm default.  The
  // particle-particle interaction table builds its Ewald coefficient from the arbitrated cutoff
  // and the direct sum tolerance (see the PPITable constructor), and the reciprocal-space
  // convolution must use the same value.  The coefficient stored in the PPPMControls object
  // follows the default &pppm cutoff, so it is correct only for a &dynamics cutoff of 8.0.
  section(1);
  const std::vector<double> cutoffs = { 8.0, 9.0, 10.0 };
  const std::vector<double> ew_coeff_ans = { 0.348639, 0.307676, 0.275106 };
  for (size_t i = 0; i < cutoffs.size(); i++) {
    const std::string cut_str = realToString(cutoffs[i], 4, 1, NumberFormat::STANDARD_REAL);
    const TextFile tf("&dynamics cut = " + cut_str + " &end", TextOrigin::RAM);
    int start_line = 0;
    const DynamicsControls dyncon(tf, &start_line, nullptr, ExceptionResponse::DIE);
    const PPPMControls pmecon;
    check(dyncon.getTranscript().getKeywordStatus("cut") == InputStatus::USER_SPECIFIED,
          "A cutoff of " + cut_str + " given in a &dynamics namelist is not recorded as "
          "user-specified.");
    const double elec_cut = arbitrateCutoff(dyncon, pmecon, NonbondedTheme::ELECTROSTATIC);
    check(elec_cut, RelationalOperator::EQUAL, cutoffs[i], "The electrostatic cutoff arbitrated "
          "between &dynamics and &pppm does not follow the &dynamics cutoff of " + cut_str + ".");
    const double direct_ew_coeff = ewaldCoefficient(elec_cut, pmecon.getDirectSumTolerance());
    const double recip_ew_coeff = arbitrateEwaldCoefficient(dyncon, pmecon);
    check(recip_ew_coeff, RelationalOperator::EQUAL, Approx(direct_ew_coeff).margin(tiny),
          "The Ewald coefficient for the reciprocal-space convolution does not match the "
          "direct-space coefficient for a &dynamics cutoff of " + cut_str + ".");
    check(recip_ew_coeff, RelationalOperator::EQUAL, Approx(ew_coeff_ans[i]).margin(1.0e-6),
          "The Ewald coefficient for the reciprocal-space convolution does not meet expectations "
          "for a &dynamics cutoff of " + cut_str + ".");
    if (fabs(cutoffs[i] - default_pme_cutoff) < tiny) {
      check(recip_ew_coeff, RelationalOperator::EQUAL,
            Approx(pmecon.getEwaldCoefficient()).margin(tiny), "The Ewald coefficient for the "
            "reciprocal-space convolution changed for a &dynamics cutoff equal to the &pppm "
            "default.");
    }
  }

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
