//
// $Id: deapActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file deapActionInitialization.hh
/// \brief Definition of the deapActionInitialization class

#ifndef deapActionInitialization_h
#define deapActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class deapActionInitialization : public G4VUserActionInitialization
{
  public:
    deapActionInitialization();
    virtual ~deapActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
