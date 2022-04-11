#include <cstdint>
#include <vector>
#include <string>
#include <tuple>
#include <ostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <functional>
#include <bitset>
#include <exception>
#include <gtest/gtest.h>
#include "../qrencode.h"
#include "../qrinput.h"
#include <quirc.h>
#include <quirc_internal.h>
#include <zbar.h>
#include <Image.h>
#include <Scanner.h>
#include <qrcodegen.hpp>
#include "compression.hpp"

extern "C" {
    extern int QRspec_getDataLength(int version, QRecLevel level);
    extern int QRspec_lengthIndicator(QRencodeMode mode, int version);
}

using check_func_t = std::function<void(const std::vector<uint8_t>&, int32_t, int32_t)>;

class CheckFunc {
public:
    CheckFunc(check_func_t func, std::string name) :
    m_name(name),
    m_func(func)
    {
    }
    ~CheckFunc() = default;

    check_func_t m_func;
    std::string m_name;
    void check (const std::vector<uint8_t>& frame, int32_t width, int32_t high) {
        m_func(frame, width, high);
    };
};

class EncodeTest : public testing::TestWithParam<std::tuple<QRecLevel, int32_t, CheckFunc>> {
public:
        EncodeTest() = default;
        virtual ~EncodeTest() = default;
protected:
    int calculateDataSize(int32_t version, QRecLevel level){
        auto data_size = QRspec_getDataLength(version, level);
        //some magic to match data size with this https://www.qrcode.com/en/about/version.html
        auto l = QRspec_lengthIndicator(QR_MODE_8, version);
        auto m = 1 << l;
        auto num = (data_size + m - 1) / m;
        auto indicator_bits = num * (MODE_INDICATOR_SIZE + l);
        indicator_bits += 7;
        data_size -= indicator_bits / 8;
        return data_size;
    }

    void generateData(int32_t version, QRecLevel level){
        std::uniform_int_distribution<uint8_t> distrib(0);
        auto data_size = calculateDataSize(version, level);
        m_random_data.resize(data_size);
        std::generate(m_random_data.begin(), m_random_data.end(), [&]{
            return distrib(m_gen);
        });
    }

    std::mt19937 m_gen{47};
    std::vector<uint8_t> m_random_data;
};

static void decodeByQuirc(const std::vector<uint8_t>& data, int32_t width, int32_t height){
    struct quirc* qr = quirc_new();
    quirc_resize(qr, width, height);
    uint8_t* p_image = quirc_begin(qr, &width, &height);
    std::copy_n(data.data(), width * height, p_image);
    quirc_end(qr);
    int32_t num_codes = quirc_count(qr);
    EXPECT_NE(0, num_codes);
    //EXPECT_EQ(1, num_codes);
    quirc_destroy(qr);
}

static void decodeByZbar(const std::vector<uint8_t>& frame, int32_t width, int32_t height){
        zbar::ImageScanner scanner;
        scanner.set_config(zbar::ZBAR_NONE, zbar::ZBAR_CFG_ENABLE, 0);
        scanner.set_config(zbar::ZBAR_QRCODE, zbar::ZBAR_CFG_ENABLE, 1);
        zbar::Image image(width, height, std::string("GREY"), NULL, width * height);
        image.set_data((void*)frame.data(), frame.size());
        int32_t decodeResult = scanner.scan(image);
        ASSERT_GT(decodeResult, 0);
}

static void fillFrame(std::vector<uint8_t>& frame, int32_t frame_width, int32_t frame_height, int32_t x_offset, int32_t y_offset,
                        uint8_t* p_qr_data ,int32_t qr_width, int32_t draw_color, int32_t bg_color, int32_t qr_scale){
    
    ASSERT_TRUE(frame_width >= x_offset + qr_width * qr_scale);
    ASSERT_TRUE(frame_height >= y_offset + qr_width * qr_scale);
    frame.clear();
    frame.resize(frame_width*frame_height, bg_color);
    auto frameIt = frame.begin() + y_offset * frame_width;

    for(int32_t y = 0; y < qr_width; y++){
        for(int32_t x = 0; x < qr_width; x++){
            uint8_t val = (0x01 & p_qr_data[x + y * qr_width]) ? draw_color : bg_color;
            std::fill_n(frameIt + x_offset + x * qr_scale, qr_scale, val);
        }
        frameIt += frame_width;
        for(uint32_t cnt = 1; cnt < qr_scale; ++cnt){
            std::copy_n(frameIt - frame_width, frame_width, frameIt);
            frameIt += frame_width;
        }
    }
}

static void dumpFrame(const std::vector<uint8_t>& frame, std::string file_name){
    std::ofstream ofs;
    ofs.open(file_name.c_str(), std::ofstream::out | std::ofstream::binary);
    ofs.write(reinterpret_cast<const char*>(frame.data()), frame.size());
    ofs.flush();
    ofs.close();
}

TEST_P(EncodeTest, LIBQRENCODE) {
    const auto& param = GetParam();
    const auto level = std::get<0>(param);
    const auto version = std::get<1>(param);
    auto data_size = calculateDataSize(version, level);

    for(int32_t i = 0; i< 100; ++i ) {
        generateData(version, level);
        auto* qr_code = QRcode_encodeData(data_size, m_random_data.data(), version, level);
        ASSERT_NE(qr_code, nullptr);
        std::vector<uint8_t> frame;
        const int32_t scale = 4;
        const int32_t fr_width = 1080;//qr_code->width*scale;
        const int32_t fr_height = 720;//qr_code->width*scale;
        const int32_t x_offset = (fr_width - qr_code->width*scale)/2;
        const int32_t y_offset = (fr_height - qr_code->width*scale)/2;
        fillFrame(frame, fr_width, fr_height, x_offset, y_offset, qr_code->data, qr_code->width, 0, 255, scale);
        dumpFrame(frame, std::to_string(version) + "_" + std::to_string(level) + ".raw" );
        ASSERT_NE(qr_code, nullptr);
        std::get<2>(param).m_func(frame, fr_width, fr_height);
    }
};

using namespace qrcodegen;

static void fillFrameQRgen(std::vector<uint8_t>& frame, int32_t frame_width, int32_t frame_height, int32_t x_offset, int32_t y_offset,
                        const QrCode& qr_code, int32_t draw_color, int32_t bg_color, int32_t qr_scale){
    
    const int32_t qrSize = qr_code.getSize();
    ASSERT_TRUE(frame_width >= x_offset + qrSize * qr_scale);
    ASSERT_TRUE(frame_height >= y_offset + qrSize * qr_scale);
    frame.clear();
    frame.resize(frame_width*frame_height, bg_color);
    auto frameIt = frame.begin() + y_offset * frame_width;

    for(int32_t y = 0; y < qrSize; y++){
        for(int32_t x = 0; x < qrSize; x++){
            uint8_t val = qr_code.getModule(x, y) ? draw_color : bg_color;
            std::fill_n(frameIt + x_offset + x * qr_scale, qr_scale, val);
        }
        frameIt += frame_width;
        for(uint32_t cnt = 1; cnt < qr_scale; ++cnt){
            std::copy_n(frameIt - frame_width, frame_width, frameIt);
            frameIt += frame_width;
        }
    }
}

static QrCode::Ecc getQRgenLevel(QRecLevel level){
    switch(level){
        case QR_ECLEVEL_L:
            return QrCode::Ecc::LOW;
        case QR_ECLEVEL_M:
            return QrCode::Ecc::MEDIUM;
        case QR_ECLEVEL_Q:
            return QrCode::Ecc::QUARTILE;
        case QR_ECLEVEL_H:
            return QrCode::Ecc::HIGH;
    }
}

TEST_P(EncodeTest, QRCODEGEN) {
    const auto& param = GetParam();
    const auto level = std::get<0>(param);
    const auto version = std::get<1>(param);
    auto data_size = calculateDataSize(version, level);

    for(int32_t i = 0; i< 100; ++i ){
        generateData(version, level);
        
        QrCode qr = QrCode::encodeBinary(m_random_data, getQRgenLevel(level));
        ASSERT_EQ(qr.getVersion(), version);
        ASSERT_EQ(qr.getErrorCorrectionLevel(), getQRgenLevel(level));
        std::vector<uint8_t> frame;
        const int32_t scale = 4;
        const int32_t fr_width = 1080;//qr_code->width*scale;
        const int32_t fr_height = 720;//qr_code->width*scale;
        const int32_t x_offset = (fr_width - qr.getSize()*scale)/2;
        const int32_t y_offset = (fr_height - qr.getSize()*scale)/2;
        fillFrameQRgen(frame, fr_width, fr_height, x_offset, y_offset, qr, 0, 255, scale);
        std::get<2>(param).m_func(frame, fr_width, fr_height);
    }
};

INSTANTIATE_TEST_SUITE_P(All, EncodeTest,
                            testing::Combine(
                                testing::Values(QR_ECLEVEL_L,
                                QR_ECLEVEL_M,
                                QR_ECLEVEL_Q,
                                QR_ECLEVEL_H
                                ),
                                testing::Range(1, QRSPEC_VERSION_MAX + 1, 1),
                                testing::ValuesIn({
                                    CheckFunc{check_func_t(decodeByQuirc), "quirc"},
                                    //CheckFunc{check_func_t(decodeByZbar), "zbar"}
                                    })
                                ),
                                [](const testing::TestParamInfo<EncodeTest::ParamType>& info){
                                    auto name = "level_" + std::to_string(std::get<0>(info.param))
                                     + "_version_" + std::to_string(std::get<1>(info.param))
                                     + "_" + std::get<2>(info.param).m_name;
                                    return name;
                                }
                         );

static void fillFrameRestoredQR(std::vector<uint8_t>& frame, int32_t frame_width, int32_t frame_height, int32_t x_offset, int32_t y_offset,
                        const std::vector<std::vector<bool>>& restored_qr, int32_t draw_color, int32_t bg_color, int32_t qr_scale){
    
    const int32_t qr_size = static_cast<int32_t>(restored_qr.size());
    ASSERT_EQ(qr_size, restored_qr.size());
    ASSERT_TRUE(frame_width >= x_offset + qr_size * qr_scale);
    ASSERT_TRUE(frame_height >= y_offset + qr_size * qr_scale);
    frame.clear();
    frame.resize(frame_width*frame_height, bg_color);
    auto frameIt = frame.begin() + y_offset * frame_width;

    for(int32_t y = 0; y < qr_size; y++){
        for(int32_t x = 0; x < qr_size; x++){
            uint8_t val = restored_qr.at(y).at(x) ? draw_color : bg_color;
            std::fill_n(frameIt + x_offset + x * qr_scale, qr_scale, val);
        }
        frameIt += frame_width;
        for(uint32_t cnt = 1; cnt < qr_scale; ++cnt){
            std::copy_n(frameIt - frame_width, frame_width, frameIt);
            frameIt += frame_width;
        }
    }
}

TEST_P(EncodeTest, QRCODEGEN_INTERCOMPRESSION) {
    const auto& param = GetParam();
    const auto level = std::get<0>(param);
    const auto version = std::get<1>(param);
    auto data_size = calculateDataSize(version, level);

    for(int32_t i = 0; i< 1; ++i ){
        generateData(version, level);
        //const auto data1 = m_random_data;
        const QrCode qr_ref = QrCode::encodeBinary(m_random_data, getQRgenLevel(level));
        generateData(version, level);
        const auto ref_data = m_random_data;
        const QrCode qr = QrCode::encodeBinary(m_random_data, getQRgenLevel(level));

        QRCompression compression;
        compression.compress(qr_ref, qr);
        std::cout << "Error: " << compression.getCompressionError() << std::endl;
        ASSERT_EQ(qr_ref.getVersion(), version);
        ASSERT_EQ(qr_ref.getErrorCorrectionLevel(), getQRgenLevel(level));
        std::vector<std::vector<bool>> restored_qr;
        compression.restoreQR(qr_ref, restored_qr);
        const auto code_size = qr_ref.getSize();
        ASSERT_EQ(restored_qr.size(), qr_ref.getSize());

        std::vector<uint8_t> frame;
        const int32_t scale = 4;
        const int32_t fr_width = 1080;//qr_code->width*scale;
        const int32_t fr_height = 720;//qr_code->width*scale;
        const int32_t x_offset = (fr_width - code_size*scale)/2;
        const int32_t y_offset = (fr_height - code_size*scale)/2;
        uint32_t err = 0;
        uint32_t counter = 0;
        for(auto y = 0; y < code_size; ++y) {
            for(auto x = 0; x < code_size; ++x) {
                ++counter;
                err += (qr.getModule(x, y) == restored_qr.at(y).at(x)) ? 0 : 1;
            }
        }

        fillFrameRestoredQR(frame, fr_width, fr_height, x_offset, y_offset, restored_qr, 0, 255, scale);
        std::get<2>(param).m_func(frame, fr_width, fr_height);
        fillFrameQRgen(frame, fr_width, fr_height, x_offset, y_offset, qr_ref, 0, 255, scale);
        std::get<2>(param).m_func(frame, fr_width, fr_height);
    }
};