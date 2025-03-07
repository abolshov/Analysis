#ifndef IMG_TFMR_HPP
#define IMG_TFMR_HPP

#include "Storage.hpp"

#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

inline constexpr Float_t HALF_LEN = 6.3;
inline constexpr Float_t PIXEL_SZ = 0.04;
inline constexpr size_t CNT = static_cast<size_t>(2*HALF_LEN/PIXEL_SZ);
inline constexpr size_t NUM_PIX = CNT % 2 == 0 ? CNT + 1 : CNT;
inline constexpr size_t ZERO_PIX = (NUM_PIX - 1)/2;

class ImageTransformer 
{
    public:
    ImageTransformer(TString output_name);
    void Transform(Storage const& storage, TTree* tree, int evt);

    private:
    // img[phi_idx][eta_idx]
    Float_t m_img_en[NUM_PIX][NUM_PIX];
    Float_t m_img_px[NUM_PIX][NUM_PIX];
    Float_t m_img_py[NUM_PIX][NUM_PIX];
    Float_t m_img_pz[NUM_PIX][NUM_PIX];
    std::unique_ptr<TFile> m_output;
    std::unique_ptr<TTree> m_tree;

    void ResetImgArr(Float_t img[NUM_PIX][NUM_PIX]);
    void DepictP4(LorentzVectorF_t const& p4, Float_t radius = 0.0);
    std::pair<size_t, size_t> FindPixelIdx(Float_t eta, Float_t phi);
    std::pair<Float_t, Float_t> FindPixelPos(size_t i, size_t j);
    std::unique_ptr<TH2F> ToHist(Float_t img[NUM_PIX][NUM_PIX], TString const& title);
};

#endif