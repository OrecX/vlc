/*****************************************************************************
 * ZrleFrameMaker.hpp
 *****************************************************************************
 * Copyright (C) 2019 OrecX LLC
 *
 * Author: Budi Kurniawan <bkurniawan@orecx.com>
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston MA 02110-1301, USA.
 *****************************************************************************/
#ifndef VLC_FBS_ZRLEFRAMEMAKER_H_
#define VLC_FBS_ZRLEFRAMEMAKER_H_
#include "zlib.h"
#include <memory>

namespace fbs {
#define MAX_INFLATE_SIZE_ZLIB 128000

struct FbsPixelFormat {
    uint8_t bitsPerPixel;
    uint8_t depth;
    uint8_t bigEndianFlag;
    uint8_t trueColorFlag;
    uint16_t redMax;
    uint16_t greenMax;
    uint16_t blueMax;
    uint8_t redShift;
    uint8_t greenShift;
    uint8_t blueShift;
};

class ZrleFrameMaker {
public:
    uint16_t frameBufferWidth;
    uint16_t frameBufferHeight;
    uint8_t bytesPerPixel;
    uint8_t zrleTilePixels24[3 * 64 * 64]; // A tile of 64 x 64 pixels, 3 bytes per pixel
    ZrleFrameMaker(demux_t *p_demux, uint16_t fbWidth, uint16_t fbHeight,
            FbsPixelFormat* fbsPixelFormat);
    void handleFrame(demux_t *p_demux, const uint8_t *data,
            const std::shared_ptr<uint8_t[]>& canvasPtr);
    void reset();
    z_stream zStream;
    unsigned long inflated;

private:
    ssize_t inf(const uint8_t* input, uint64_t inputLength, uint8_t *buf);
    void populateColorArray(uint8_t* data, uint64_t* pos, uint8_t colorArray[],
            uint16_t paletteSize);
    void handleRawPixels(uint8_t* data, uint64_t* pos, uint8_t tileWidth,
            uint8_t tileHeight);
    void handlePackedPixels(uint8_t* data, uint64_t* pos, uint8_t tileWidth,
            uint8_t tileHeight, uint8_t colorArray[], uint16_t paletteSize);
    void handlePlainRLEPixels(uint8_t* data, uint64_t* pos, uint8_t tileWidth,
            uint8_t tileHeight);
    void handlePackedRLEPixels(uint8_t* data, uint64_t* pos, uint8_t tileWidth,
            uint8_t tileHeight, uint8_t colorArray[]);
    void copyTileToBuffer(uint8_t *canvas, uint32_t j, uint32_t i,
            uint8_t tileWidth, uint8_t tileHeight);
    uint32_t readPixel(uint8_t* data, uint64_t* pos);
};
} // namespace
#endif /* VLC_FBS_ZRLEFRAMEMAKER_H_ */
