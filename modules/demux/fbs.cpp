/*****************************************************************************
 * fbs.cpp : fbs demuxer
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
#include <memory>
#include <vlc_common.h>
#include <vlc_demux.h>
#include <vlc_fs.h>
#include <vlc_plugin.h>
#include "zlib.h"

#define MAX_INFLATE_SIZE_ZLIB 128000
/** Initial memory alignment of data block.
 * @note This must be a multiple of sizeof(void*) and a power of two.
 * libavcodec AVX optimizations require at least 32-bytes. */
#define BLOCK_ALIGN        32

/** Initial reserved header and footer size. */
#define BLOCK_PADDING      32
/*****************************************************************************
 * Module descriptor
 *****************************************************************************/
static int Open(vlc_object_t *);
static void Close(vlc_object_t *);
static int Seek(demux_t *, vlc_tick_t);
vlc_module_begin ()
    set_shortname( "FBS" )
    set_description( N_("FBS stream demuxer" ) )
    set_capability( "demux", 212 )
    set_callbacks( Open, Close )
    set_category( CAT_INPUT )
    set_subcategory( SUBCAT_INPUT_DEMUX )
vlc_module_end ()

/*****************************************************************************
 * Definitions of structures used by this plugin
 *****************************************************************************/
typedef struct {
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
} FbsPixelFormat;

typedef struct {
    uint64_t frame_size;
    es_out_id_t *p_es_video;
    es_format_t fmt_video;
    date_t pcr;
    uint8_t rfbVersion;
    uint16_t frameBufferWidth;
    uint16_t frameBufferHeight;
    uint8_t framesPerSecond;
    uint16_t headerPos;
    uint32_t timestamp;
    uint32_t lastTimestamp;
    uint64_t canvasLength;
    bool seeking;
    FbsPixelFormat fbsPixelFormat;
    uint64_t inflated;
    z_stream zStream;
    uint8_t zrleTilePixels24[3 * 64 * 64]; // A tile of 64 x 64 pixels, 3 bytes per pixel
    uint8_t *canvas;
} demux_sys_t;

static int Demux(demux_t *);
static int Control(demux_t *, int, va_list);

void copyTileToBuffer(const demux_sys_t *p_sys, const uint8_t bitsPerPixel, const uint16_t frameBufferWidth,
        const uint32_t j, const uint32_t i, const uint8_t tileWidth, const uint8_t tileHeight);
void handleFrame(demux_sys_t *p_sys, const uint8_t bitsPerPixel,
        const uint16_t frameBufferWidth, const uint8_t *data);
ssize_t inf(demux_sys_t *p_sys, const uint8_t* input, const uint64_t inputLength, uint8_t *buf);
void handlePackedPixels(uint8_t *zrleTilePixels24, const uint8_t* data, uint64_t* pos,
        const uint8_t tileWidth, const uint8_t tileHeight, const uint8_t colorArray[], const uint16_t paletteSize);
void handlePackedRLEPixels(uint8_t *zrleTilePixels24, const uint8_t* data, uint64_t* pos, const uint8_t tileWidth,
        const uint8_t tileHeight, const uint8_t colorArray[]);
void handlePlainRLEPixels(uint8_t *zrleTilePixels24, const uint8_t bitsPerPixel, const uint8_t* data,
        uint64_t* pos, const uint8_t tileWidth, const uint8_t tileHeight);
void populateColorArray(const uint8_t bitsPerPixel, const uint8_t* data, uint64_t* pos,
        uint8_t colorArray[], const uint16_t paletteSize);
int readLastTimestamp(const char *);
uint32_t readPixel(const uint8_t bitsPerPixel, const uint8_t* data, uint64_t* pos);
void resetZStream(demux_sys_t *p_sys);
void updateCanvas(const demux_t *p_demux, const uint8_t bitsPerPixel,
        const uint16_t frameBufferWidth, const uint64_t frameSize, const uint32_t);

/*****************************************************************************
 * Open: initializes FBS demux structures
 *****************************************************************************/
static int Open(vlc_object_t *p_this) {
    demux_t *p_demux = (demux_t*) p_this;
    demux_sys_t *p_sys;
    p_demux->p_sys = p_sys = (demux_sys_t*) malloc(sizeof(demux_sys_t));
    const vlc_fourcc_t i_chroma = VLC_CODEC_RGB24;
    const uint8_t i_sar_num = 1;
    const uint8_t i_sar_den = 1;
    const uint8_t *p_peek;

    p_sys->lastTimestamp = readLastTimestamp(p_demux->psz_filepath);
    p_sys->framesPerSecond = 5;
    p_sys->timestamp = 0;
    p_sys->seeking = false;
    int readBytes = vlc_stream_Peek(p_demux->s, &p_peek, 2000);
    if (readBytes == -1) {
        return VLC_EGENERIC;
    }
    if (strncmp((char*) p_peek, "FBS 001.000", 11)) {
        return VLC_EGENERIC; // file invalid
    }

    uint64_t pos = 12;
    for (uint8_t frameNo = 1; frameNo < 6; frameNo++) {
        int dataLength = U32_AT(&p_peek[pos]);
        int paddedDataLength = 4 * ((dataLength + 3) / 4);
        if (frameNo == 1) {
            p_sys->rfbVersion = p_peek[pos + 14] - 48; // rfbVersion is either 3 or 8
        }
        if (frameNo == 3 && p_sys->rfbVersion == 3) {
            // no security result for RFB 3.3, skip
            continue;
        }
        if (frameNo == 4) {
            p_sys->frameBufferWidth = U16_AT(&p_peek[pos + 4]);
            p_sys->frameBufferHeight = U16_AT(&p_peek[pos + 6]);
            p_sys->fbsPixelFormat.bitsPerPixel = p_peek[pos + 8];
            p_sys->fbsPixelFormat.depth = p_peek[pos + 9];
            p_sys->fbsPixelFormat.bigEndianFlag = p_peek[pos + 10];
            p_sys->fbsPixelFormat.trueColorFlag = p_peek[pos + 11];
            p_sys->fbsPixelFormat.redMax = U16_AT(&p_peek[pos + 12]);
            p_sys->fbsPixelFormat.greenMax = U16_AT(&p_peek[pos + 14]);
            p_sys->fbsPixelFormat.blueMax = U16_AT(&p_peek[pos + 16]);
            p_sys->fbsPixelFormat.redShift = p_peek[pos + 18];
            p_sys->fbsPixelFormat.greenShift = p_peek[pos + 19];
            p_sys->fbsPixelFormat.blueShift = p_peek[pos + 20];
        }
        pos += 4 + paddedDataLength + /* timestamp */4;
        p_sys->headerPos = pos; // used by SEEK
    }

    /* Set the demux function */
    es_format_Init(&p_sys->fmt_video, VIDEO_ES, i_chroma);
    video_format_Setup(&p_sys->fmt_video.video, i_chroma,
            p_sys->frameBufferWidth, p_sys->frameBufferHeight,
            p_sys->frameBufferWidth, p_sys->frameBufferHeight,
            i_sar_num, i_sar_den);

    date_Init(&p_sys->pcr, p_sys->fmt_video.video.i_frame_rate,
            p_sys->fmt_video.video.i_frame_rate_base);
    date_Set(&p_sys->pcr, VLC_TICK_0);

    const vlc_chroma_description_t *dsc = vlc_fourcc_GetChromaDescription(
            p_sys->fmt_video.video.i_chroma);
    p_sys->frame_size = 0; // need to be set to 0
    for (unsigned i = 0; i < dsc->plane_count; i++) {
        unsigned pitch = (p_sys->frameBufferWidth + (dsc->p[i].w.den - 1))
                * dsc->p[i].w.num / dsc->p[i].w.den * dsc->pixel_size;
        unsigned lines = (p_sys->frameBufferHeight + (dsc->p[i].h.den - 1))
                * dsc->p[i].h.num / dsc->p[i].h.den;
        p_sys->frame_size += pitch * lines;
    }
    p_sys->p_es_video = es_out_Add(p_demux->out, &p_sys->fmt_video);
    p_demux->pf_demux = Demux;
    p_demux->pf_control = Control;
    p_sys->canvasLength = p_sys->frameBufferWidth * p_sys->frameBufferHeight * 3;
    resetZStream(p_sys);
    // skip to data
    readBytes = vlc_stream_Read(p_demux->s, NULL, pos);
    p_sys->canvas = new uint8_t[p_sys->canvasLength];
    return VLC_SUCCESS;
}

/*****************************************************************************
 * Close: frees unused data
 *****************************************************************************/
static void Close(vlc_object_t *p_this) {
    demux_t *p_demux = (demux_t*) p_this;
    demux_sys_t *p_sys = (demux_sys_t*) p_demux->p_sys;
    delete[] p_sys->canvas;
    delete p_sys;
}

/*****************************************************************************
 * Control:
 *****************************************************************************/
static int Control(demux_t *p_demux, int i_query, va_list args) {
    demux_sys_t *p_sys = (demux_sys_t*) p_demux->p_sys;
    const int64_t i_bps = 8LL * p_sys->canvasLength;
    switch (i_query) {
    case DEMUX_GET_LENGTH:
        //assign a number in microseconds: timestamp (in ms) * 1000
        *va_arg(args, vlc_tick_t *) = p_sys->lastTimestamp * 1000;
        return VLC_SUCCESS;
    case DEMUX_CAN_SEEK:
        *va_arg(args, bool *) = true;
        return VLC_SUCCESS;
    case DEMUX_GET_TIME:
        *va_arg(args, vlc_tick_t *) = p_sys->timestamp * 1000;
        return VLC_SUCCESS;
    case DEMUX_GET_POSITION:
        double *pf;
        pf = va_arg(args, double *);
        *pf = (double) p_sys->timestamp / (double) p_sys->lastTimestamp;
        return VLC_SUCCESS;
    case DEMUX_SET_POSITION:
        double f = va_arg(args, double);
        vlc_tick_t i64 = f * 1000;
        return Seek(p_demux, i64);
    }
    return demux_vaControlHelper(p_demux->s, 0, -1, i_bps,
            p_sys->canvasLength, i_query, args);
}

/*****************************************************************************
 * Demux: reads and demuxes data packets
 *****************************************************************************
 * Returns -1 in case of error, 0 in case of EOF, 1 otherwise
 *****************************************************************************/
static int Demux(demux_t *p_demux) {
    demux_sys_t *p_sys = (demux_sys_t *) p_demux->p_sys;
    vlc_tick_t pcr = date_Get(&p_sys->pcr);
    es_out_SetPCR(p_demux->out, pcr);
    if (!p_sys->seeking) {
        updateCanvas(p_demux, p_sys->fbsPixelFormat.bitsPerPixel,
                p_sys->frameBufferWidth, p_sys->frame_size, pcr / 1000);
    }
    p_sys->pcr.i_divider_num = p_sys->framesPerSecond; //how many times in a second Demux() is called
    p_sys->pcr.i_divider_den = 1;
    date_Increment(&p_sys->pcr, 1);
    return VLC_DEMUXER_SUCCESS;
}

static int Seek(demux_t *p_demux, vlc_tick_t i_date) {
    demux_sys_t *p_sys = (demux_sys_t *) p_demux->p_sys;
    uint32_t soughtTimestamp = i_date * p_sys->lastTimestamp / 1000;
    // Rewind to right after the header
    if (vlc_stream_Seek(p_demux->s, p_sys->headerPos) == -1) {
        return VLC_DEMUXER_EOF;
    }
    p_sys->seeking = true;
    p_sys->pcr.i_divider_num = 0;
    resetZStream(p_sys);

    p_sys->timestamp = 0;
    updateCanvas(p_demux, p_sys->fbsPixelFormat.bitsPerPixel,
            p_sys->frameBufferWidth, p_sys->frame_size, soughtTimestamp);
    p_sys->seeking = false;
    p_sys->pcr.date = soughtTimestamp * 1000;
    return VLC_SUCCESS;
}

void copyTileToBuffer(const demux_sys_t *p_sys, const uint8_t bitsPerPixel, const uint16_t frameBufferWidth,
        const uint32_t x, const uint32_t y, const uint8_t tileWidth, const uint8_t tileHeight) {
    assert(tileWidth <= 64 && tileHeight <= 64);
    if (bitsPerPixel == 32) {
        for (uint8_t i = 0; i < tileHeight; i++) {
            for (uint8_t j = 0; j < tileWidth; j++) {
                int source = 3 * (i * tileWidth + j);
                int position = ((y + i) * frameBufferWidth + (x + j)) * 3;
                p_sys->canvas[position + 2] = p_sys->zrleTilePixels24[source];     //red;
                p_sys->canvas[position + 1] = p_sys->zrleTilePixels24[source + 1]; //green;
                p_sys->canvas[position] = p_sys->zrleTilePixels24[source + 2];     //blue;
            }
        }
    }
}

void handleFrame(demux_sys_t *p_sys, const uint8_t bitsPerPixel,
        const uint16_t frameBufferWidth, const uint8_t *data) {
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
        const uint32_t rectDataLength = U32_AT(&data[pos]);
        pos += 4;

        uint8_t* inflatedData = new uint8_t[10 * MAX_INFLATE_SIZE_ZLIB];
        if (!inflatedData) {
            return;
        }
        inf(p_sys, data + pos, rectDataLength, inflatedData);

        pos += rectDataLength;
        uint64_t inflatedDataReadPosition = 0;
        for (int j = rectY; j < rectY + rectHeight; j += 64) {
            // last row tiles may have height that is < 64
            int tileHeight = std::min(64, rectY + rectHeight - j);
            for (int k = rectX; k < rectX + rectWidth; k += 64) {
                // last column tiles may have width that is < 64
                int tileWidth = std::min(64, rectX + rectWidth - k);
                uint8_t subencoding = inflatedData[inflatedDataReadPosition++];
                bool runLength = (subencoding & 0x80) != 0;
                uint8_t paletteSize = subencoding & 0x7F;
                uint8_t colorArray[128 * 3];
                populateColorArray(bitsPerPixel, inflatedData,
                        &inflatedDataReadPosition, colorArray, paletteSize);
                // read palette colors
                if (paletteSize == 1) {
                    for (uint16_t d = j; d < j + tileHeight; d++) {
                        for (uint16_t e = k; e < k + tileWidth; e++) {
                            uint32_t index = (d * frameBufferWidth + e) * 3;
                            p_sys->canvas[index++] = colorArray[0]; //red;
                            p_sys->canvas[index++] = colorArray[1]; //green;
                            p_sys->canvas[index] = colorArray[2];   //blue;
                        }
                    }
                    continue;
                }
                if (!runLength) {
                    if (paletteSize == 0) {
                        populateColorArray(bitsPerPixel, inflatedData, &inflatedDataReadPosition,
                                p_sys->zrleTilePixels24, tileWidth * tileHeight);
                    } else {
                        handlePackedPixels(p_sys->zrleTilePixels24, inflatedData, &inflatedDataReadPosition,
                                tileWidth, tileHeight, colorArray, paletteSize);
                    }
                } else {
                    if (paletteSize == 0) {
                        handlePlainRLEPixels(p_sys->zrleTilePixels24, bitsPerPixel, inflatedData,
                                &inflatedDataReadPosition, tileWidth, tileHeight);
                    } else {
                        handlePackedRLEPixels(p_sys->zrleTilePixels24, inflatedData, &inflatedDataReadPosition, tileWidth,
                                tileHeight, colorArray);
                    }
                }
                copyTileToBuffer(p_sys, bitsPerPixel, frameBufferWidth, k, j, tileWidth, tileHeight);
            }
        }
        delete[] inflatedData;
    }
}

void handlePackedPixels(uint8_t *zrleTilePixels24, const uint8_t* data, uint64_t* pos, const uint8_t tileWidth,
        const uint8_t tileHeight, const uint8_t colorArray[], const uint16_t paletteSize) {
    assert(tileWidth <= 64 && tileHeight <= 64);
    uint8_t shift = 1;
    if (paletteSize > 16) {
        shift = 8;
    } else if (paletteSize > 4) {
        shift = 4;
    } else if (paletteSize > 2) {
        shift = 2;
    }
    uint16_t index = 0;
    for (uint8_t i = 0; i < tileHeight; ++i) {
        uint16_t end = index + tileWidth;
        uint8_t counter1;
        uint16_t counter2 = 0;
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

void handlePackedRLEPixels(uint8_t *zrleTilePixels24, const uint8_t* data, uint64_t* pos,
        const uint8_t tileWidth, const uint8_t tileHeight, const uint8_t colorArray[]) {
    assert(tileWidth <= 64 && tileHeight <= 64);
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

void handlePlainRLEPixels(uint8_t *zrleTilePixels24, const uint8_t bitsPerPixel,
        const uint8_t* data, uint64_t* pos, const uint8_t tileWidth, const uint8_t tileHeight) {
    assert(tileWidth <= 64 && tileHeight <= 64);
    uint16_t index = 0;
    uint16_t end = tileWidth * tileHeight;
    while (index < end) {
        uint8_t length;
        uint32_t pixel = readPixel(bitsPerPixel, data, pos);
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

ssize_t inf(demux_sys_t *p_sys, const uint8_t* input, const uint64_t inputLength, uint8_t *buf) {
    p_sys->zStream.avail_in = inputLength;
    p_sys->zStream.next_in = (uint8_t*) input;
    uint8_t chunkData[MAX_INFLATE_SIZE_ZLIB];
    uint32_t count = 0;
    while(1) {
        p_sys->zStream.avail_out = MAX_INFLATE_SIZE_ZLIB;
        p_sys->zStream.next_out = chunkData;
        int inflateResult = inflate(&p_sys->zStream, Z_NO_FLUSH);
        uint32_t numBytesLeft = p_sys->zStream.total_out - p_sys->inflated;
        for (unsigned int i = 0; i < numBytesLeft; i++) {
            buf[count++] = chunkData[i];
        }
        p_sys->inflated = p_sys->zStream.total_out;
        if (inflateResult != Z_OK || p_sys->zStream.avail_in == 0
                || (p_sys->zStream.avail_out != 0 && inflateResult == Z_STREAM_END)) {
            break;
        }
    }
    return count;
}

void populateColorArray(const uint8_t bitsPerPixel, const uint8_t* rectData, uint64_t* pos,
        uint8_t colorArray[], const uint16_t paletteSize) {
    assert(paletteSize <= 64 * 64);
    if (bitsPerPixel == 32) {
        for (uint16_t i = 0; i < paletteSize; i++) {
            uint16_t index = i * 3;
            colorArray[index++] = rectData[(*pos)++]; // red
            colorArray[index++] = rectData[(*pos)++]; // green
            colorArray[index] = rectData[(*pos)++]; // blue
        }
    }
}

int readLastTimestamp(const char* filepath) {
    unsigned char c[4]; // read last 4 bytes
    FILE *file = vlc_fopen(filepath, "rb");
    if (!file) {
        printf("Error opening file\n");
        return -1;
    }
    fseek(file, -4, SEEK_END);
    int bytesRead = fread(c, sizeof(char), 4, file);
    fclose(file);
    if (bytesRead == 4) {
        return U32_AT(&c);
    } else {
        return -1;
    }
}

uint32_t readPixel(const uint8_t bitsPerPixel, const uint8_t* data, uint64_t* pos) {
    uint32_t pixel = 0;
    if (bitsPerPixel == 32) {
        uint8_t red = data[(*pos)++];
        uint8_t green = data[(*pos)++];
        uint8_t blue = data[(*pos)++];
        pixel = red << 16 | green << 8 | blue;
    }
    return pixel;
}

void resetZStream(demux_sys_t *p_sys) {
    p_sys->inflated = p_sys->zStream.avail_in = 0;
    p_sys->zStream.zalloc = Z_NULL;
    p_sys->zStream.zfree = Z_NULL;
    p_sys->zStream.opaque = Z_NULL;
    p_sys->zStream.next_in = Z_NULL;
    inflateInit(&p_sys->zStream);
}

void updateCanvas(const demux_t *p_demux, const uint8_t bitsPerPixel,
        const uint16_t frameBufferWidth, const uint64_t frameSize, const uint32_t t) {
    demux_sys_t *p_sys = (demux_sys_t *) p_demux->p_sys;
    while (p_sys->timestamp <= t) {
        block_t *p_block = vlc_stream_Block(p_demux->s, 4);
        if (p_block == NULL) {
            return;
        }
        int dataLength = U32_AT(&p_block->p_buffer[0]);
        int paddedDataLength = 4 * ((dataLength + 3) / 4);
        block_Release(p_block);
        p_block = vlc_stream_Block(p_demux->s, paddedDataLength + /*timestamp*/4);
        p_block->i_size = frameSize + BLOCK_ALIGN + 2 * BLOCK_PADDING;
        p_block->i_buffer = frameSize;
        p_sys->timestamp = U32_AT(&p_block->p_buffer[paddedDataLength]);
        handleFrame(p_sys, bitsPerPixel, frameBufferWidth, p_block->p_buffer);
        block_Release(p_block);
    }
    block_t *p_block = block_Alloc(p_sys->frame_size);
    p_block->i_buffer = p_sys->frame_size;
    memcpy(p_block->p_buffer, p_sys->canvas, p_sys->frame_size);
    p_block->i_dts = p_block->i_pts = p_sys->timestamp * 1000;
    es_out_Send(p_demux->out, p_sys->p_es_video, p_block);
}
