#include "ImageTransformer.hpp"

#include <stdexcept>

ImageTransformer::ImageTransformer(TString output_name)
:   m_output(std::make_unique<TFile>(output_name, "RECREATE"))
,   m_tree(std::make_unique<TTree>("img_tree", "img_tree"))
{
    ResetImgArr(m_img_en);
    ResetImgArr(m_img_px);
    ResetImgArr(m_img_py);
    ResetImgArr(m_img_pz);
    m_tree->Branch("image_en", m_img_en, "m_img_en[NUM_PIX][NUM_PIX]/F");
    m_tree->Branch("image_px", m_img_px, "m_img_px[NUM_PIX][NUM_PIX]/F");
    m_tree->Branch("image_py", m_img_py, "m_img_py[NUM_PIX][NUM_PIX]/F");
    m_tree->Branch("image_pz", m_img_pz, "m_img_pz[NUM_PIX][NUM_PIX]/F");
}

void ImageTransformer::ResetImgArr(Float_t img[NUM_PIX][NUM_PIX])
{
    for (size_t r = 0; r < NUM_PIX; ++r)
    {
        for (size_t c = 0; c < NUM_PIX; ++c)
        {
            img[r][c] = 0.0;
        }
    }
}

std::pair<size_t, size_t> ImageTransformer::FindPixelIdx(Float_t eta, Float_t phi)
{
    // phi is periodic, eta is not
    // x axis (eta) is indexing columns
    // y axis (phi) is indexin rows
    size_t col = static_cast<size_t>((eta + HALF_LEN)/PIXEL_SZ);
    size_t row = static_cast<size_t>((phi + HALF_LEN)/PIXEL_SZ) % NUM_PIX;
    return {row, col};
}

std::pair<Float_t, Float_t> ImageTransformer::FindPixelPos(size_t i, size_t j)
{
    Float_t y = -1.0*HALF_LEN + PIXEL_SZ/2.0 + PIXEL_SZ*i;
    Float_t x = -1.0*HALF_LEN + PIXEL_SZ/2.0 + PIXEL_SZ*j;
    return {x, y};
}

void ImageTransformer::DepictP4(LorentzVectorF_t const& p4, Float_t radius)
{
    Float_t phi = p4.Phi();
    Float_t eta = p4.Eta();
    Float_t r2 = radius*radius;

    auto [row, col] = FindPixelIdx(phi, eta);
    
    if (radius > 0.0)
    {
        // draw object of non-zero radius
        Float_t min_phi = phi - radius;
        Float_t min_eta = eta - radius;

        // at first need to count all pixels that need to be painted
        size_t n = static_cast<size_t>(radius/PIXEL_SZ);
        size_t num_steps = 2*n + 1;
        int filled_pixel_count = 0;
        for (size_t i = 0; i < num_steps; ++i)
        {
            Float_t cur_eta = min_eta + i*PIXEL_SZ;
            if (std::abs(cur_eta) > HALF_LEN)
            {
                continue;
            }

            for (size_t j = 0; j < num_steps; ++j)
            {
                Float_t cur_phi = min_phi + j*PIXEL_SZ;
                Float_t dsit2 = (cur_phi - phi)*(cur_phi - phi) + (cur_eta - eta)*(cur_eta - eta);
                if (dsit2 <= r2)
                {
                    // count this pixel
                    ++filled_pixel_count;
                }
            }
        }

        // paint pixels in separate loop
        // don't know how to optimize it yet
        for (size_t i = 0; i < num_steps; ++i)
        {
            Float_t cur_eta = min_eta + i*PIXEL_SZ;
            if (std::abs(cur_eta) > HALF_LEN)
            {
                continue;
            }

            for (size_t j = 0; j < num_steps; ++j)
            {
                Float_t cur_phi = min_phi + j*PIXEL_SZ;
                Float_t dsit2 = (cur_phi - phi)*(cur_phi - phi) + (cur_eta - eta)*(cur_eta - eta);
                if (dsit2 <= r2)
                {
                    auto [r, c] = FindPixelIdx(cur_eta, cur_phi); // this function needs to be fixed; now it does UB
                    m_img_en[r][c] += p4.E()/filled_pixel_count;
                    m_img_px[r][c] += p4.Px()/filled_pixel_count;
                    m_img_py[r][c] += p4.Py()/filled_pixel_count;
                    m_img_pz[r][c] += p4.Pz()/filled_pixel_count;
                }
            }
        }
    }
    else 
    {
        // draw object of radius size (i.e. 1 pixel)
        m_img_en[row][col] += p4.E();
        m_img_px[row][col] += p4.Px();
        m_img_py[row][col] += p4.Py();
        m_img_pz[row][col] += p4.Pz();
    }
} 

void ImageTransformer::Transform(Storage const& storage, TTree* tree, int evt)
{
    tree->GetEntry(evt);
    // TODO
    m_tree->Fill();
}

std::unique_ptr<TH2F> ImageTransformer::ToHist(Float_t img[NUM_PIX][NUM_PIX], TString const& title)
{
    auto hist = std::make_unique<TH2F>("img", title, NUM_PIX, -HALF_LEN, HALF_LEN, NUM_PIX, -HALF_LEN, HALF_LEN);
    hist->SetStats(0);
    for (size_t r = 0; r < NUM_PIX; ++r)
    {
        for (size_t c = 0; c < NUM_PIX; ++c)
        {
            int bin = hist->FindBin(img[r][c]);
            hist->SetBinContent(bin, img[r][c]);
        }
    }
    return hist;
}