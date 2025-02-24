#ifndef EST_RCDR_HPP
#define EST_RCDR_HPP

#include <memory>

#include "TFile.h"
#include "TTree.h"

class EstimationRecorder
{
    public:
    EstimationRecorder(TString const& dbg_file_name);

    inline void ResetTree(std::unique_ptr<TTree> tree) { m_tree = std::move(tree); }
    inline void ReleaseTree() { m_tree.release(); }
    inline void FillTree() { m_tree->Fill(); }
    inline void WriteTree(TString const& tree_name) { m_file->WriteObject(m_tree.get(), tree_name); }

    private:
    std::unique_ptr<TFile> m_file;
    std::unique_ptr<TTree> m_tree;
};

#endif