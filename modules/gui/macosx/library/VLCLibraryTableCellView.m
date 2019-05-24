/*****************************************************************************
 * VLCLibraryTableCellView.m: MacOS X interface module
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

#import "VLCLibraryTableCellView.h"
#import "extensions/NSFont+VLCAdditions.h"
#import "views/VLCImageView.h"
#import "views/VLCTrackingView.h"
#import "main/VLCMain.h"
#import "library/VLCLibraryController.h"
#import "library/VLCLibraryDataTypes.h"

@interface VLCLibraryTableCellView ()
{
    VLCLibraryController *_libraryController;
}
@end

@implementation VLCLibraryTableCellView

- (void)awakeFromNib
{
    self.singlePrimaryTitleTextField.font = [NSFont VLClibraryCellTitleFont];
    self.primaryTitleTextField.font = [NSFont VLClibraryCellTitleFont];
    self.secondaryTitleTextField.font = [NSFont VLClibraryCellSubtitleFont];
    [self prepareForReuse];
}

- (void)prepareForReuse
{
    [super prepareForReuse];
    self.representedImageView.image = nil;
    self.primaryTitleTextField.hidden = YES;
    self.secondaryTitleTextField.hidden = YES;
    self.singlePrimaryTitleTextField.hidden = YES;
    self.trackingView.viewToHide = nil;
    self.playInstantlyButton.hidden = YES;
}

- (void)setRepresentedMediaItem:(VLCMediaLibraryMediaItem *)representedMediaItem
{
    _representedMediaItem = representedMediaItem;
    self.trackingView.viewToHide = self.playInstantlyButton;
}

- (IBAction)playInstantly:(id)sender
{
    if (!_libraryController) {
        _libraryController = [[VLCMain sharedInstance] libraryController];
    }

    [_libraryController appendItemToPlaylist:_representedMediaItem playImmediately:YES];
}

@end
