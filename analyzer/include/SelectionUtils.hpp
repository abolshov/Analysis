#ifndef SELEC_UTILS_HPP
#define SELEC_UTILS_HPP

#ifdef DEV
    #include "Event.hpp"
    #include "Constants.hpp"

    // true if event is passing a selection, false otherwise
    bool IsRecoverable(Event const& s, Channel ch);
    bool IsFiducial(Event const& s, VecLVF_t const& jets, Channel ch);

    bool CorrRecoLep(int lep_type, int lep_genLep_kind);
#endif 

#endif