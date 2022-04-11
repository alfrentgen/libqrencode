#ifndef QRCOMPRESSION_H
#define QRCOMPRESSION_H

#include <vector>
#include <tuple>
#include <functional>
#include <qrcodegen.hpp>

constexpr uint8_t MAX_QR_SIZE = 177;

/*enum class path_t {
    LINE,
    MEANDER,
};

enum class direction_t {
    FORWARD,
    BACKWARD,
};

enum class comparison_t {
    //comp_to_ref
    ROW_TO_ROW,
    COL_TO_ROW,
    ROW_TO_COL,
    COL_TO_COL,
};*/

typedef struct {
    uint8_t x;
    uint8_t y;
} coordinate_t;

typedef struct {
    int8_t  x_shift : 3;
    int8_t  y_shift : 3;
    uint8_t step : 2;
} path_t;

static constexpr int8_t MAX_SHIFT = 3;
static constexpr int8_t MIN_SHIFT = -3;
static constexpr uint8_t MAX_STEP = 3;

typedef struct {
    coordinate_t coord;
    path_t path;
} ref_line_descriptor_t;

using namespace qrcodegen;
using line_compression_t = std::tuple<uint8_t, ref_line_descriptor_t, uint64_t>;
using compress_func_t = std::function<void(const QrCode& ref, const QrCode& comp, std::vector<line_compression_t>& line_compressions)>;

class QRCompression {
public:
    QRCompression() = default;
    ~QRCompression() = default;

    void compress(const QrCode& ref, const QrCode& comp);
    void restoreQR(const QrCode& ref, std::vector<std::vector<bool>>& restored_qr);
    uint64_t getCompressionError() const;
    const std::vector<line_compression_t>& getLineCompressions() const {
        return m_compressions;
    };
private:
    std::vector<line_compression_t> m_compressions;
    static std::map<uint8_t, compress_func_t> m_compress_functions;
    
    compress_func_t& getCompressionFunction(uint8_t qr_size);
};

#endif