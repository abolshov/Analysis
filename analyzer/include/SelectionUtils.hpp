#ifndef SELEC_UTILS_HPP
#define SELEC_UTILS_HPP

#include "Storage.hpp"

// true if event is passing a selection, false otherwise
bool IsRecoverable(Storage const& s, Channel ch, bool top_sel = true);
bool IsFiducial(Storage const& s, VecLVF_t const& jets, Channel ch);

bool CorrRecoLep(int lep_type, int lep_genLep_kind);

#endif