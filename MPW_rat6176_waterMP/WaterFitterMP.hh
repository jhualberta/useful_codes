////////////////////////////////////////////////////////////////////////
/// \class RAT::WaterFitter
///
/// \brief   Water phase fitter
///
/// \author  Phil Jones <p.jones22@physics.ox.ac.uk>
/// \author Matt Mottram < m.mottram@qmul.ac.uk> -- contact person
///
/// REVISION HISTORY:
///     - 27/04/2011 : P.Jones - First Revision, new file.
///     - 2014-03-29 : P G Jones - Updated lifetime, added BeginOfRun method
///     - 2014-07-20 : A Soerenson - Added muon fitting.
///     - 2014-08-27 : M Mottram - Updated fit components as per DocDB 2157
///     - 2014-12-09 : P G Jones - Allow chained PMTSelectors
///     - 2015-02-24 : E Beier - added Isotropy classifier, see DocDB 2835 \n
///     - 2015-03-23 : M Mottram - Fix incase no FitVertex in seed
///     - 2015-04-08 : M Mottram - removed Beta14
///     - 2015-05-29 : M Mottram - use energyPromptLookup method
///     - 2016-01-14 : M Mottram - Add SetFOMs
///     - 2016-01-18 : J Walker - Change energy method to energyRSP
///     - 2016-03-15 : M Mottram - Fail waterFitter if quad fails
///     - 2017-05-22 : J Caravaca - Add PMTCal selector
///
/// \details  Best (as defined by reconstruction group) water phase fitter
///          combination.
///          Will FAIL if either position/time or energy components
///          throw.
///
////////////////////////////////////////////////////////////////////////

#ifndef __RAT_WaterFitterMP__
#define __RAT_WaterFitterMP__

#include <RAT/Processor.hh>

namespace RAT
{
namespace DS
{
  class EV;
}
namespace Methods
{
  class Method;
}
namespace PDFs
{
  class PDF;
}
namespace Optimisers
{
  class Optimiser;
}
namespace PMTSelectors
{
  class PMTSelector;
}
namespace Classifiers
{
  class Classifier;
}
namespace DS{
  class FitResult;
}

class WaterFitterMP : public Processor
{
public:
  WaterFitterMP();
  virtual ~WaterFitterMP();

  /// Called at the beginning of runs to setup the waterFitter
  ///
  /// @param[in] run
  void BeginOfRun( DS::Run& run );

  /// Called to invoke the waterFitter
  ///
  /// @param[in, out] run Data structure for the event
  /// @param[in, out] ds Data structure for the event
  virtual Processor::Result DSEvent( DS::Run& run, DS::Entry& ds );

  /// Called at the end of runs
  ///
  /// @param[in] run
  void EndOfRun( DS::Run& run );
protected:

  virtual Processor::Result Event( DS::Run& run, DS::EV& ev );

  Processor::Result SetResult( DS::FitResult& waterResult, DS::EV& ev, size_t currentPass);
  void DriveCorrection( DS::FitResult& fitResult );

  // Default fitters may have subroutines that set the same FOM name for different parameters
  // Add a prefix for the FOM name in the final result to distinguish them.
  void SetFOMs( DS::FitResult& fitResult, const DS::FitResult& partialResult, const std::string& prefix );

  Methods::Method* fQuadSeed;
  Methods::Method*  fPositionTime;
  Optimisers::Optimiser* fMetaDrivePowell;
  PMTSelectors::PMTSelector* fModeCut;
  PDFs::PDF* fGV1D;

  Methods::Method* fEnergySeed;
  Methods::Method* fEnergyMethod;

  Methods::Method* fMuonWater;

  Methods::Method* fDirectionSeed;
  Methods::Method* fDirection;
  Methods::Method*  fPositionTimeDirection;
  Optimisers::Optimiser* fSimulatedAnnealing;
  PMTSelectors::PMTSelector* fDirectionModeCut;
  PMTSelectors::PMTSelector* fTimeResidualCut0;
  PMTSelectors::PMTSelector* fTimeResidualCut;
  PMTSelectors::PMTSelector* fPMTCalSelector;
  PDFs::PDF* fPositionDirectionPDF;

  Classifiers::Classifier* fIsotropy;
  Classifiers::Classifier* fITR;
  Classifiers::Classifier* fQPDT;

  Methods::Method* fSmearResult;
  void SetI(const std::string& param, const int value);
  void SetD(const std::string& param, const double value);
  void SetS(const std::string& param, const std::string& value);

  size_t fQuadNhitCut; ///< Maximum nhit for which quad method cannot run
  double fDrivePositionScale;
  double fDriveDirectionScale;
  int fDriveEnable;
};

} // ::RAT

#endif
