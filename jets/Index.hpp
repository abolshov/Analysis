#ifndef INDEX_HPP
#define INDEX_HPP

#include <iostream>
#include "TROOT.h"

// struct to save indices of the bbWW decay chain
struct GenPartIndex
{
    Int_t h1 = -1;
    Int_t h2 = -1;
    Int_t w1 = -1;
    Int_t w2 = -1;
    Int_t b1 = -1;
    Int_t b2 = -1;
    Int_t q1 = -1;
    Int_t q2 = -1;
    Int_t l = -1;
    Int_t nu = -1;

    inline operator bool() const noexcept
    {
        if (h1 != -1 && h2 != -1 && w1 != -1 && w2 != -1 && b1 != -1 && b2 != -1 && q1 != -1 && q2 != -1 && l != -1 && nu != -1) return true;
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, GenPartIndex const& idx);
};

struct GenJetIndex
{
    Int_t bj1 = -1;
    Int_t bj2 = -1;
    Int_t lj1 = -1;
    Int_t lj2 = -1;

    inline operator bool() const noexcept
    {
        if (bj1 != -1 && bj2 != -1 && bj1 != bj2 && lj1 != -1 && lj2 != -1 && lj1 != lj2) return true;
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, GenJetIndex const& idx);
};

struct PtEtaPhiMArray
{
    Float_t* ptPtr = nullptr;
    Float_t* etaPtr = nullptr;
    Float_t* phiPtr = nullptr;
    Float_t* mPtr = nullptr;
    UInt_t n = 0;

    inline operator bool() const noexcept 
    {
        return (ptPtr && etaPtr && phiPtr && mPtr);
    } 
};

#endif