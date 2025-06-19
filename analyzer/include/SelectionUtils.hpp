#ifndef SELEC_UTILS_HPP
#define SELEC_UTILS_HPP

#ifdef DEV
    #include "Event.hpp"
    #include "Constants.hpp"

    // true if event is passing a selection, false otherwise
    bool IsRecoverable(Event const& s, Channel ch, Topology bb_top, Topology qq_top = Topology::Resolved);
    bool IsFiducial(Event const& s, VecLVF_t const& jets, Channel ch, Topology bb_top, Topology qq_top = Topology::Resolved);

    bool CorrRecoLep(int lep_type, int lep_genLep_kind);
#endif 

#endif