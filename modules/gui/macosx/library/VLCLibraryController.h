/*****************************************************************************
 * VLCLibraryController.h: MacOS X interface module
 *****************************************************************************
 * Copyright (C) 2019 VLC authors and VideoLAN
 *
 * Authors: Felix Paul Kühne <fkuehne # videolan -dot- org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA 02110-1301, USA.
 *****************************************************************************/

#import <Foundation/Foundation.h>

@class VLCLibraryModel;
@class VLCMediaLibraryMediaItem;

NS_ASSUME_NONNULL_BEGIN

@interface VLCLibraryController : NSObject

@property (readonly) VLCLibraryModel *libraryModel;

- (int)appendItemToPlaylist:(VLCMediaLibraryMediaItem *)mediaItem playImmediately:(BOOL)playImmediately;
- (void)showItemInFinder:(VLCMediaLibraryMediaItem *)mediaItem;
- (int)attemptToGenerateThumbnailForMediaItem:(VLCMediaLibraryMediaItem *)mediaItem;

- (int)addFolderWithFileURL:(NSURL *)fileURL;
- (int)banFolderWithFileURL:(NSURL *)fileURL;
- (int)unbanFolderWithFileURL:(NSURL *)fileURL;
- (int)removeFolderWithFileURL:(NSURL *)fileURL;

- (int)clearHistory;

@end

NS_ASSUME_NONNULL_END
