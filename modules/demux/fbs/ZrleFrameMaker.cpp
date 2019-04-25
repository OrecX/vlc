/*****************************************************************************
 * ZrleFrameMaker.cpp
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
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <assert.h>
#include "vlc_common.h"
#include "ZrleFrameMaker.hpp"

namespace fbs {
ZrleFrameMaker::ZrleFrameMaker(demux_t *p_demux, uint16_t fbWidth, uint16_t fbHeight,
        FbsPixelFormat *fbsPixelFormat) :
                frameBufferWidth(fbWidth), frameBufferHeight(fbHeight),
                bytesPerPixel(fbsPixelFormat->bitsPerPixel / 8) {
    reset();
}

void ZrleFrameMaker::handleFrame(demux_t *p_demux, const uint8_t *data,
        const std::shared_ptr<uint8_t[]>& canvasPtr) {
    // byte 1 = message type, byte 2 = padding, both can be ignored
    uint16_t rectCount = U16_AT(&data[2]);
    uint64_t pos = 4;
    for (uint16_t i = 0; i < rectCount; i++) {
        uint16_t rectX = U16_AT(&data[pos]);
        pos += 2;
        uint16_t rectY = U16_AT(&data[pos]);
        pos += 2;
        uint16_t rectWidth = U16_AT(&data[pos]);
        pos += 2;
        uint16_t rectHeight = U16_AT(&data[pos]);
        pos += 2;
        pos += 4; // skip encoding
        uint32_t rectDataLength = U32_AT(&data[pos]);
        pos += 4;

        uint8_t* inflatedData = (uint8_t*) malloc(10 * MAX_INFLATE_SIZE_ZLIB);
        if (!inflatedData) {
            return;
        }
        ssize_t size = inf(data + pos, rectDataLength, inflatedData);

        pos += rectDataLength;
        uint64_t inflatedDataReadPosition = 0;
        uint8_t *canvas = canvasPtr.get();
        for (int j = rectY; j < rectY + rectHeight; j += 64) {
            // last row tiles may have height that is < 64
            int tileHeight = std::min(64, rectY + rectHeight - j);
            for (int k = rectX; k < rectX + rectWidth; k += 64) {
                // last column tiles may have width that is < 64
                int tileWidth = std::min(64, rectX + rectWidth - k);
                // the first byte of the data contains subencoding
                uint8_t subencoding = inflatedData[inflatedDataReadPosition++];
                bool runLength = (subencoding & 0x80) != 0;
                uint8_t paletteSize = subencoding & 0x7F;
                uint8_t colorArray[128 * 3];
                populateColorArray(inflatedData, &inflatedDataReadPosition, colorArray, paletteSize);
                // read palette colors
                if (paletteSize == 1) {
                    for (int d = j; d < j + tileHeight; d++) {
                        for (int e = k; e < k + tileWidth; e++) {
                            int index = (d * frameBufferWidth + e) * 3;
                            canvas[index++] = colorArray[0]; //red;
                            canvas[index++] = colorArray[1]; //green;
                            canvas[index] = colorArray[2];   //blue;
                        }
                    }
                    continue;
                }
                if (!runLength) {
                    if (paletteSize == 0) {
                        handleRawPixels(inflatedData, &inflatedDataReadPosition, tileWidth, tileHeight);
                    } else {
                        handlePackedPixels(inflatedData, &inflatedDataReadPosition, tileWidth,
                                tileHeight, colorArray, paletteSize);
                    }
                } else {
                    if (paletteSize == 0) {
                        handlePlainRLEPixels(inflatedData, &inflatedDataReadPosition, tileWidth, tileHeight);
                    } else {
                        handlePackedRLEPixels(inflatedData, &inflatedDataReadPosition, tileWidth,
                                tileHeight, colorArray);
                    }
                }
                copyTileToBuffer(canvas, k, j, tileWidth, tileHeight);
            }
        }
        delete inflatedData;
    }
}

void ZrleFrameMaker::reset() {
    inflated = zStream.avail_in = 0;
    zStream.zalloc = Z_NULL;
    zStream.zfree = Z_NULL;
    zStream.opaque = Z_NULL;
    zStream.next_in = Z_NULL;
    inflateInit(&zStream);
}

void ZrleFrameMaker::populateColorArray(uint8_t* rectData, uint64_t* pos,
        uint8_t colorArray[], uint16_t paletteSize) {
    assert(paletteSize <= 64 * 64);
    switch (bytesPerPixel) {
    case 1: // not supported
        break;
    case 2: // not supported
        break;
    case 4:
        for (uint16_t i = 0; i < paletteSize; i++) {
            uint16_t index = i * 3;
            colorArray[index++] = rectData[(*pos)++]; // red
            colorArray[index++] = rectData[(*pos)++]; // green
            colorArray[index] = rectData[(*pos)++]; // blue
        }
        break;
    }
}

void ZrleFrameMaker::handleRawPixels(uint8_t* rectData, uint64_t* pos,
        uint8_t tileWidth, uint8_t tileHeight) {
    populateColorArray(rectData, pos, zrleTilePixels24,
            tileWidth * tileHeight);
}

void ZrleFrameMaker::handlePackedPixels(uint8_t* data, uint64_t* pos,
        uint8_t tileWidth, uint8_t tileHeight, uint8_t colorArray[], uint16_t paletteSize) {
    uint8_t shift = 1;
    if (paletteSize > 16) {
        shift = 8;
    } else if (paletteSize > 4) {
        shift = 4;
    } else if (paletteSize > 2) {
        shift = 2;
    }
    int index = 0;
    for (int i = 0; i < tileHeight; ++i) {
        int end = index + tileWidth;
        uint8_t counter1;
        int counter2 = 0;
        while (index < end) {
            if (counter2 == 0) {
                counter1 = data[(*pos)++];
                counter2 = 8;
            }
            counter2 -= shift;
            uint16_t colorIndex = 3 * ((counter1 >> counter2) & ((1 << shift) - 1));
            zrleTilePixels24[index * 3] = colorArray[colorIndex++];
            zrleTilePixels24[index * 3 + 1] = colorArray[colorIndex++];
            zrleTilePixels24[index * 3 + 2] = colorArray[colorIndex];
            index++;
        }
    }
}

void ZrleFrameMaker::handlePlainRLEPixels(uint8_t* data, uint64_t* pos,
        uint8_t tileWidth, uint8_t tileHeight) {
    assert(tileWidth <= 64 && tileHeight <= 64);
    uint16_t index = 0;
    uint16_t end = tileWidth * tileHeight;
    while (index < end) {
        uint8_t length;
        uint32_t pixel = readPixel(data, pos);
        uint16_t totalLength = 1;
        do {
            length = data[(*pos)++];
            totalLength += length;
        } while (length == 255);
        assert(totalLength <= end - index);
        while (totalLength-- > 0) {
            zrleTilePixels24[3 * index] = (pixel & 0xFF0000) >> 16;
            zrleTilePixels24[3 * index + 1] = (pixel & 0xFF00) >> 8;
            zrleTilePixels24[3 * index + 2] = pixel & 0xFF;
            index++;
        }
    }
}

void ZrleFrameMaker::handlePackedRLEPixels(uint8_t* data, uint64_t* pos,
        uint8_t tileWidth, uint8_t tileHeight, uint8_t colorArray[]) {
    uint16_t index = 0;
    uint16_t end = tileWidth * tileHeight;
    while (index < end) {
        uint8_t flag = data[(*pos)++];
        uint16_t totalLength = 1;
        if ((flag & 128) != 0) {
            uint8_t length;
            do {
                length = data[(*pos)++];
                totalLength += length;
            } while (length == 255);
            assert(totalLength <= end - index);
        }
        flag &= 127;
        while (totalLength-- != 0) {
            zrleTilePixels24[3 * index] = colorArray[3 * flag];
            zrleTilePixels24[3 * index + 1] = colorArray[3 * flag + 1];
            zrleTilePixels24[3 * index + 2] = colorArray[3 * flag + 2];
            index++;
        }
    }
}

void ZrleFrameMaker::copyTileToBuffer(uint8_t *canvas,
        uint32_t x, uint32_t y, uint8_t tileWidth, uint8_t tileHeight) {
    switch (bytesPerPixel) {
    case 1: // not supported
        break;
    case 2: // not supported
        break;
    case 4:
        for (int i = 0; i < tileHeight; i++) {
            for (int j = 0; j < tileWidth; j++) {
                int source = 3 * (i * tileWidth + j);
                int position = ((y + i) * frameBufferWidth + (x + j)) * 3;
                canvas[position + 2] = zrleTilePixels24[source];     //red;
                canvas[position + 1] = zrleTilePixels24[source + 1]; //green;
                canvas[position] = zrleTilePixels24[source + 2];     //blue;
            }
        }
    }
}

uint32_t ZrleFrameMaker::readPixel(uint8_t* data, uint64_t* pos) {
    uint32_t pixel = 0;
    switch (bytesPerPixel) {
    case 1: // not supported
        break;
    case 2: // not supported
        break;
    case 4:
        uint8_t red = data[(*pos)++];
        uint8_t green = data[(*pos)++];
        uint8_t blue = data[(*pos)++];
        pixel = red << 16 | green << 8 | blue;
        break;
    }
    return pixel;
}

ssize_t ZrleFrameMaker::inf(const uint8_t* input, uint64_t inputLength, uint8_t *buf) {
    zStream.avail_in = inputLength;
    zStream.next_in = (uint8_t*) input;
    uint8_t chunkData[MAX_INFLATE_SIZE_ZLIB];
    uint32_t count = 0;
    while(1) {
        zStream.avail_out = MAX_INFLATE_SIZE_ZLIB;
        zStream.next_out = chunkData;
        int inflateResult = inflate(&zStream, Z_NO_FLUSH);
        unsigned int numBytesLeft = zStream.total_out - inflated;
        for (unsigned int i = 0; i < numBytesLeft; i++) {
            buf[count++] = chunkData[i];
        }
        inflated = zStream.total_out;
        if (inflateResult != Z_OK || zStream.avail_in == 0
                || (zStream.avail_out != 0 && inflateResult == Z_STREAM_END)) {
            break;
        }
    }
    return count;
}
} // namespace
