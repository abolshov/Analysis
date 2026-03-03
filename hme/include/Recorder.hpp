#ifndef RCDR_HPP
#define RCDR_HPP

#include <memory>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

class Recorder
{
    public:
    Recorder() = default;
    explicit Recorder(TString const& file_name);

    inline void ResetTree(std::unique_ptr<TTree> tree) { m_tree = std::move(tree); }
    inline void FillTree() { m_tree->Fill(); }
    inline void WriteTree(TString const& tree_name) { m_file->WriteObject(m_tree.get(), tree_name); }
    inline bool ShouldRecord() { return m_file != nullptr; }
    inline std::unique_ptr<TTree>& GetTree() { return m_tree; }
    void OpenFile(TString const& file_name);

    private:
    std::unique_ptr<TFile> m_file;
    std::unique_ptr<TTree> m_tree;
};

#endif