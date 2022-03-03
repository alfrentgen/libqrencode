#include <cstdint>
#include <vector>
#include <tuple>
#include <ostream>
#include <random>
#include <algorithm>
#include <gtest/gtest.h>
#include "../qrencode.h"
#include "../qrinput.h"

extern "C" {
    extern int QRspec_getDataLength(int version, QRecLevel level);
    extern int QRspec_lengthIndicator(QRencodeMode mode, int version);
}

class EncodeTest : public testing::TestWithParam<std::tuple<QRecLevel, int32_t>> {
public:
        EncodeTest() = default;
        virtual ~EncodeTest() = default;
protected:
    int calculateDataSize(int32_t version, QRecLevel level){
        auto data_size = QRspec_getDataLength(version, level);
        //some magic to match results with this https://www.qrcode.com/en/about/version.html
        auto l = QRspec_lengthIndicator(QR_MODE_8, version);
        auto m = 1 << l;
        auto num = (data_size + m - 1) / m;
        auto indicator_bits = num * (MODE_INDICATOR_SIZE + l);
        indicator_bits += 7;
        data_size -= indicator_bits / 8;
        return data_size;
    }

    std::vector<uint8_t> m_random_data;
};

TEST_P(EncodeTest, generateAll) {
    const auto& param = GetParam();
    const auto level = std::get<0>(param);
    const auto version = std::get<1>(param);
    auto data_size = calculateDataSize(version, level);

    std::cout << level << ", " << version << ": " << data_size << std::endl;
    auto* qr_input = QRinput_new2(version, QR_ECLEVEL_H);

    std::mt19937 gen(47);
    std::uniform_int_distribution<uint8_t> distrib(0);
    m_random_data.resize(data_size);
    std::generate(m_random_data.begin(), m_random_data.end(), [&]{
        return distrib(gen);
    });
    auto* qr_code = QRcode_encodeData(data_size, m_random_data.data(), version, level);
    ASSERT_NE(qr_code, nullptr);
};

INSTANTIATE_TEST_SUITE_P(All, EncodeTest, testing::Combine(
                                            testing::Values(QR_ECLEVEL_L, 
                                            QR_ECLEVEL_M, 
                                            QR_ECLEVEL_Q, 
                                            QR_ECLEVEL_H),
                                            testing::Range(1, QRSPEC_VERSION_MAX + 1, 1))
                         );
