#include <cstdint>
#include <vector>
#include <string>
#include <tuple>
#include <ostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <exception>
#include <cassert>
#include <map>
#include <utility>
#include "../qrencode.h"
#include "../qrinput.h"
#include <gtest/gtest.h>
#include "compression.hpp"

template<typename T>
class CyclicCounter {
public:
    CyclicCounter(T a, T b, T val = 0) :
    m_min {a < b ? a : b},
    m_max {b > a ? b : a},
    m_val {val}
    {
    };

    ~CyclicCounter() = default;

    T operator()() const
    {
        return m_val;
    };

    T operator++()
    {
        if(m_max == m_val++) {
            m_val = m_min;
        }
        return m_val;
    };

    T operator+=(T inc)
    {
        while(inc != 0){
            inc > 0 ? operator++(), --inc : operator--(), ++inc;
        }
        return m_val;
    };

    T operator--()
    {
        if(m_min == m_val--) {
            m_val = m_max;
        }
        return m_val;
    };

    T operator-=(T dec)
    {
        while(dec != 0) {
            dec > 0 ? operator--(), --dec : operator++(), ++dec;
        }
        return m_val;
    };

private:
    T m_min;
    T m_max;
    T m_val;
};

static void getPath(const coordinate_t start_coord, const path_t path_desc, std::vector<coordinate_t>& coords, uint8_t qr_size)
{
    ASSERT_TRUE(start_coord.x <= qr_size && start_coord.y <= qr_size);
    coords.clear();
    CyclicCounter<uint8_t> x_coord(0, qr_size-1, start_coord.x);
    CyclicCounter<uint8_t> y_coord(0, qr_size-1, start_coord.y);

    const auto step = path_desc.step + 1;
    const uint16_t path_length = qr_size * step;
    const auto x_shift_abs = std::abs(path_desc.x_shift);
    const auto y_shift_abs = std::abs(path_desc.y_shift);
    ASSERT_TRUE(x_shift_abs != 0 || y_shift_abs != 0);
    uint16_t counter{0};
    while(coords.size() < qr_size) {
        for(uint8_t x = 0; x < x_shift_abs && coords.size() != qr_size; ++x){
            if(counter % step == 0) {
                coords.push_back({x_coord(), y_coord()});
            }
            path_desc.x_shift > 0 ? ++x_coord : --x_coord;
            ++counter;
        }
        for(uint8_t y = 0; y < y_shift_abs && coords.size() != qr_size; ++y){
            if(counter % step == 0) {
                coords.push_back({x_coord(), y_coord()});
            }
            path_desc.y_shift > 0 ? ++y_coord : --y_coord;
            ++counter;
        }
    }
}

template<uint8_t QR_SIZE>
static void getRefLine(const QrCode& ref_qr, const std::vector<coordinate_t> coords, std::bitset<QR_SIZE>& ref_line)
{
    ASSERT_LE(QR_SIZE, MAX_QR_SIZE);
    ASSERT_EQ(ref_line.size(), QR_SIZE);
    ASSERT_EQ(ref_line.size(), coords.size());
    for(uint8_t i {0}; i < coords.size(); ++i) {
        ref_line[i] = ref_qr.getModule(coords[i].x, coords[i].y);
    }
}

using ECCLevel = qrcodegen::QrCode::Ecc;
static double getECCThreshold(const ECCLevel& ecc) {
    switch(ecc){
        case ECCLevel::HIGH:
        return 0.3;
        case ECCLevel::QUARTILE:
        return 0.25;
        case ECCLevel::MEDIUM:
        return 0.15;
        case ECCLevel::LOW:
        return 0.07;
    }
    EXPECT_TRUE(false);
    throw std::invalid_argument("Wrong ECC level value");
}

static void printCoords(const std::vector<coordinate_t> coords){
    for(auto& c : coords){
        std::cout << "("<< std::to_string(c.x) << ", " << std::to_string(c.y) << "); ";
    }
    std::cout << std::endl;
}

template<uint8_t QR_SIZE>
static void compress(const QrCode& ref, const QrCode& comp, std::vector<line_compression_t>& line_compressions)
{
    ASSERT_EQ(ref.getSize(), comp.getSize());
    ASSERT_EQ(ref.getSize(), QR_SIZE);
    const auto error_threshold = getECCThreshold(comp.getErrorCorrectionLevel());
    for(uint8_t comp_y {0}; comp_y < QR_SIZE; ++comp_y) {
        std::bitset<QR_SIZE> comp_line;
        for(uint8_t comp_x {0}; comp_x < QR_SIZE; ++comp_x) {
            comp_line[comp_x] = comp.getModule(comp_x, comp_y);
        }
        uint64_t lowest_error = std::numeric_limits<uint64_t>::max();
        ref_line_descriptor_t best_descriptor{};
        std::vector<coordinate_t> best_coords;
        std::bitset<QR_SIZE> ref_line;
        const auto stop_threshold = std::ceil(error_threshold * QR_SIZE);
        [&]() {
            for(uint8_t ref_y {0}; ref_y < QR_SIZE; ++ref_y) {
                for(uint8_t ref_x {0}; ref_x < QR_SIZE; ++ref_x) {
                    for(int8_t y_shift {MIN_SHIFT}; y_shift <= MAX_SHIFT; ++y_shift) {
                        for(int8_t x_shift {MIN_SHIFT}; x_shift <= MAX_SHIFT; ++x_shift) {
                            if(y_shift == 0 && x_shift == 0) { // path collapses into a dot, skip
                                continue;
                            }
                            for(uint8_t step {0}; step <= MAX_STEP; ++step) {
                                const ref_line_descriptor_t ref_line_desc{{ref_x, ref_y}, {x_shift, y_shift, step}};
                                std::vector<coordinate_t> coords;
                                getPath(ref_line_desc.coord, ref_line_desc.path, coords, QR_SIZE);
                                ASSERT_EQ(coords.size(), QR_SIZE);
                                getRefLine<QR_SIZE>(ref, coords, ref_line);
                                const uint64_t error = (comp_line^ref_line).count();
                                ASSERT_LE(error, QR_SIZE);
                                if(error < lowest_error) {
                                    lowest_error = error;
                                    best_descriptor = ref_line_desc;
                                    best_coords = coords;
                                    if(error < stop_threshold) {
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }();
        //printCoords(best_coords);
        line_compressions.emplace_back(comp_y, best_descriptor, lowest_error);
    }
};

std::map<uint8_t, compress_func_t> QRCompression::m_compress_functions;

template <std::uint8_t... is>
static void fillCompressFunctions(std::integer_sequence<uint8_t, is...>, std::map<uint8_t, compress_func_t>& compress_functions)
{
    ([&](){compress_functions[21 + is*4] = compress_func_t{&compress<21 + is*4>};}(), ...);
}

compress_func_t& QRCompression::getCompressionFunction(uint8_t qr_size)
{
    return QRCompression::m_compress_functions.at(qr_size);
}

void QRCompression::compress(const QrCode& ref, const QrCode& comp)
{
    if(m_compressions.size()) {
        m_compressions.clear();
    }
    if(QRCompression::m_compress_functions.empty()) {
        fillCompressFunctions(std::make_integer_sequence<uint8_t, 40>{}, m_compress_functions);
    }
    getCompressionFunction(ref.getSize())(ref, comp, m_compressions);
}

uint64_t QRCompression::getCompressionError() const {
    uint64_t error{0};
    for(const auto& comp : m_compressions){
        error += std::get<2>(comp);
    }
    return error;
}

void QRCompression::restoreQR(const QrCode& ref, std::vector<std::vector<bool>>& restored_qr){
    const auto qr_size = ref.getSize();
    ASSERT_EQ(qr_size, m_compressions.size());
    ASSERT_EQ(qr_size, m_compressions.size());
    restored_qr.resize(qr_size, std::vector<bool>(qr_size));
    std::vector<coordinate_t> coords;
    size_t counter{0};
    for(const auto& line_comp : m_compressions){
        const auto& ref_line_desc = std::get<1>(line_comp);
        getPath(ref_line_desc.coord, ref_line_desc.path, coords, qr_size);
        ASSERT_EQ(qr_size, coords.size());
        //printCoords(coords);
        for(const auto& coord : coords) {
            ASSERT_LT(coord.x, qr_size);
            ASSERT_LT(coord.y, qr_size);
            /*if (coord.x > qr_size || coord.y > qr_size) {
                printf("%d %d %d\n", coord.x, coord.y, qr_size);
            }*/
            restored_qr.at(counter/qr_size).at(counter%qr_size) = ref.getModule(coord.x, coord.y);
            ++counter;
        }
    }
}