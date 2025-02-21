#ifndef EST_RCDR_HPP
#define EST_RCDR_HPP

#include <memory>

#include "TFile.h"
#include "TTree.h"

class EstimationRecorder
{
    public:
    explicit EstimationRecorder(TString const& out_file_name);

    inline void ResetTree(std::unique_ptr<TTree> tree) { m_tree = std::move(tree); }
    inline void ReleaseTree() { m_tree.release(); }
    inline void WriteTree() { m_file->Write(); }
    inline void FillTree() { m_tree->Fill(); }

    private:
    std::unique_ptr<TFile> m_file;
    std::unique_ptr<TTree> m_tree;
};

#endif