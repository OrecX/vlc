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
#include <vlc_common.h>
#include <vlc_plugin.h>
#include <vlc_demux.h>
#include <fcntl.h>
#include "ZrleFrameMaker.hpp"

/** Initial memory alignment of data block.
 * @note This must be a multiple of sizeof(void*) and a power of two.
 * libavcodec AVX optimizations require at least 32-bytes. */
#define BLOCK_ALIGN        32

/** Initial reserved header and footer size. */
#define BLOCK_PADDING      32
/*****************************************************************************
 * Module descriptor
 *****************************************************************************/
namespace fbs {
static int Open(vlc_object_t *);
static void Close(vlc_object_t *);
static int Seek(demux_t *, vlc_tick_t);
} // namespace

using namespace fbs;
vlc_module_begin ()
set_shortname( "FBS" )
set_description( N_("FBS stream demuxer" ) )
set_capability( "demux", 212 )
set_callbacks( Open, Close )
set_category( CAT_INPUT )
set_subcategory( SUBCAT_INPUT_DEMUX )
vlc_module_end ()

namespace fbs {
/*****************************************************************************
 * Definitions of structures used by this plugin
 *****************************************************************************/
typedef struct {
    int frame_size;
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
} demux_sys_t;

static FbsPixelFormat *fbsPixelFormat;
static ZrleFrameMaker *zrleFrameMaker;
static block_t *p_block;
static std::shared_ptr<uint8_t[]> canvasPtr;
static int Demux(demux_t *);
static int Control(demux_t *, int, va_list);
static int readLastTimestamp(char *);
/*****************************************************************************
 * Open: initializes FBS demux structures
 *****************************************************************************/
static int Open(vlc_object_t * p_this) {
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
    if (strncmp((char*) &p_peek[0], "FBS 001.000", 11)) {
        // file invalid
        return VLC_EGENERIC;
    }

    std::string header((char*) &p_peek[0], 2000);

    fbsPixelFormat = new FbsPixelFormat();
    uint64_t pos = 12;
    for (int frameNo = 1; frameNo < 6; frameNo++) {
        int dataLength = U32_AT(&p_peek[pos]);
        std::string data = header.substr(pos + 4, dataLength);
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
            fbsPixelFormat->bitsPerPixel = p_peek[pos + 8];
            fbsPixelFormat->depth = p_peek[pos + 9];
            fbsPixelFormat->bigEndianFlag = p_peek[pos + 10];
            fbsPixelFormat->trueColorFlag = p_peek[pos + 11];
            fbsPixelFormat->redMax = U16_AT(&p_peek[pos + 12]);
            fbsPixelFormat->greenMax = U16_AT(&p_peek[pos + 14]);
            fbsPixelFormat->blueMax = U16_AT(&p_peek[pos + 16]);
            fbsPixelFormat->redShift = p_peek[pos + 18];
            fbsPixelFormat->greenShift = p_peek[pos + 19];
            fbsPixelFormat->blueShift = p_peek[pos + 20];
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
    canvasPtr = std::shared_ptr<uint8_t[]>(new uint8_t[p_sys->canvasLength],
            std::default_delete<uint8_t[]>());

    zrleFrameMaker = new ZrleFrameMaker(p_demux,
            p_sys->frameBufferWidth, p_sys->frameBufferHeight, fbsPixelFormat);
    // skip to data
    readBytes = vlc_stream_Read(p_demux->s, NULL, pos);
    return VLC_SUCCESS;
}

/*****************************************************************************
 * Close: frees unused data
 *****************************************************************************/
static void Close(vlc_object_t *p_this) {
    demux_t *p_demux = (demux_t*) p_this;
    demux_sys_t *p_sys = (demux_sys_t*) p_demux->p_sys;
    delete p_sys;
    delete p_block;
    delete fbsPixelFormat;
    delete zrleFrameMaker;
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
    vlc_tick_t i_pcr = date_Get(&p_sys->pcr);
    es_out_SetPCR(p_demux->out, i_pcr);

    while (p_sys->timestamp <= i_pcr / 1000 && !p_sys->seeking) {
        p_block = vlc_stream_Block(p_demux->s, 4);
        if (p_block == NULL) {
            return VLC_DEMUXER_EOF;
        }
        int dataLength = U32_AT(&p_block->p_buffer[0]);
        int paddedDataLength = 4 * ((dataLength + 3) / 4);
        p_block = vlc_stream_Block(p_demux->s,
                paddedDataLength + /*timestamp*/4);
        p_block->i_size = p_sys->frame_size + BLOCK_ALIGN + 2 * BLOCK_PADDING;
        p_block->i_buffer = p_sys->frame_size;
        p_sys->timestamp = U32_AT(&p_block->p_buffer[paddedDataLength]);
        zrleFrameMaker->handleFrame(p_demux, p_block->p_buffer, canvasPtr);
        p_block->p_buffer = canvasPtr.get();
        p_block->i_dts = p_block->i_pts = i_pcr;
        es_out_Send(p_demux->out, p_sys->p_es_video, p_block);
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
    int readBytes = vlc_stream_Seek(p_demux->s, p_sys->headerPos);
    if (readBytes == -1) {
        return VLC_DEMUXER_EOF;
    }
    p_sys->seeking = true;
    p_sys->pcr.i_divider_num = 0;
    zrleFrameMaker->reset();

    p_sys->timestamp = 0;
    while (p_sys->timestamp <= soughtTimestamp) {
        p_block = vlc_stream_Block(p_demux->s, 4);
        if (p_block == NULL) {
            return VLC_DEMUXER_EOF;
        }
        uint32_t dataLength = U32_AT(p_block->p_buffer);
        uint32_t paddedDataLength = 4 * ((dataLength + 3) / 4);
        p_block = vlc_stream_Block(p_demux->s,
                paddedDataLength + /*timestamp*/4);
        p_block->i_size = p_sys->frame_size + BLOCK_ALIGN + 2 * BLOCK_PADDING;
        p_block->i_buffer = p_sys->frame_size;
        p_sys->timestamp = U32_AT(&p_block->p_buffer[paddedDataLength]);
        zrleFrameMaker->handleFrame(p_demux, p_block->p_buffer, canvasPtr);
        p_block->p_buffer = canvasPtr.get();
        p_block->i_dts = p_block->i_pts = p_sys->timestamp * 1000; //i_pcr;
        es_out_Send(p_demux->out, p_sys->p_es_video, p_block);
    }
    p_sys->seeking = false;
    p_sys->pcr.date = soughtTimestamp * 1000;
    return VLC_SUCCESS;
}

int readLastTimestamp(char* filepath) {
    int fd;
    unsigned char c[4]; // read last 4 bytes
    fd = open(filepath, O_RDONLY);
    if (fd == -1) {
        printf("Error opening file\n");
        return -1;
    }
    lseek(fd, -4L, SEEK_END);
    int readBytes = read(fd, c, 4); // Read 4 bytes
    close(fd);
    if (readBytes == 4) {
        return U32_AT(&c);
    } else {
        return -1;
    }
}
} // namespace
