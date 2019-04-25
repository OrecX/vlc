/*****************************************************************************
 * VLCLibraryWindow.m: MacOS X interface module
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

#import "VLCLibraryWindow.h"
#import "extensions/NSString+Helpers.h"
#import "library/VLCLibraryCollectionViewItem.h"
#import "main/VLCMain.h"

#import "playlist/VLCPlaylistTableCellView.h"
#import "playlist/VLCPlaylistController.h"
#import "playlist/VLCPlaylistDataSource.h"

#import "windows/mainwindow/VLCControlsBarCommon.h"
#import "windows/video/VLCFSPanelController.h"
#import "windows/video/VLCVoutView.h"

static const float f_min_window_width = 604.;
static const float f_min_window_height = 307.;
static const float f_playlist_row_height = 72.;

static NSString *VLCLibraryCellIdentifier = @"VLCLibraryCellIdentifier";

@interface VLCLibraryWindow ()
{
    VLCPlaylistDataSource *_playlistDataSource;
    VLCLibraryDataSource *_libraryDataSource;

    NSRect _windowFrameBeforePlayback;
}
@end

@implementation VLCLibraryWindow

- (void)awakeFromNib
{
    _fspanel = [[VLCFSPanelController alloc] init];
    [_fspanel showWindow:self];

    _segmentedTitleControl.segmentCount = 3;
    [_segmentedTitleControl setTarget:self];
    [_segmentedTitleControl setAction:@selector(segmentedControlAction)];
    [_segmentedTitleControl setLabel:_NS("Music") forSegment:0];
    [_segmentedTitleControl setLabel:_NS("Video") forSegment:1];
    [_segmentedTitleControl setLabel:_NS("Network") forSegment:2];
    [_segmentedTitleControl sizeToFit];

    VLCPlaylistController *playlistController = [[VLCMain sharedInstance] playlistController];
    _playlistDataSource = [[VLCPlaylistDataSource alloc] init];
    _playlistDataSource.playlistController = playlistController;
    _playlistDataSource.tableView = _playlistTableView;
    playlistController.playlistDataSource = _playlistDataSource;

    _playlistTableView.dataSource = _playlistDataSource;
    _playlistTableView.delegate = _playlistDataSource;
    _playlistTableView.rowHeight = f_playlist_row_height;
    [_playlistTableView reloadData];

    _libraryDataSource = [[VLCLibraryDataSource alloc] init];
    _libraryCollectionView.dataSource = _libraryDataSource;
    _libraryCollectionView.delegate = _libraryDataSource;
    [_libraryCollectionView registerClass:[VLCLibraryCollectionViewItem class] forItemWithIdentifier:VLCLibraryCellIdentifier];
    [_libraryCollectionView reloadData];
}

- (void)segmentedControlAction
{
}

- (void)playlistDoubleClickAction:(id)sender
{
    NSInteger selectedRow = self.playlistTableView.selectedRow;
    if (selectedRow == -1)
        return;

    [[[VLCMain sharedInstance] playlistController] playItemAtIndex:selectedRow];
}

- (void)changePlaylistState:(int)event
{
}

- (void)videoplayWillBeStarted
{
    if (!self.fullscreen)
        _windowFrameBeforePlayback = [self frame];
}

- (void)setVideoplayEnabled
{
    BOOL b_videoPlayback = [[VLCMain sharedInstance] activeVideoPlayback];

    if (!b_videoPlayback) {
        if (!self.nonembedded && (!self.nativeFullscreenMode || (self.nativeFullscreenMode && !self.fullscreen)) && _windowFrameBeforePlayback.size.width > 0 && _windowFrameBeforePlayback.size.height > 0) {

            // only resize back to minimum view of this is still desired final state
            CGFloat f_threshold_height = f_min_video_height + [self.controlsBar height];
            if (_windowFrameBeforePlayback.size.height > f_threshold_height) {
                if ([[VLCMain sharedInstance] isTerminating])
                    [self setFrame:_windowFrameBeforePlayback display:YES];
                else
                    [[self animator] setFrame:_windowFrameBeforePlayback display:YES];

            }
        }

        _windowFrameBeforePlayback = NSMakeRect(0, 0, 0, 0);

        [self makeFirstResponder: _playlistTableView];
        [[[VLCMain sharedInstance] voutProvider] updateWindowLevelForHelperWindows: NSNormalWindowLevel];

        // restore alpha value to 1 for the case that macosx-opaqueness is set to < 1
        [self setAlphaValue:1.0];
        [self.videoView setHidden:YES];
    } else {
        [self.videoView setHidden:NO];
    }

    if (self.nativeFullscreenMode) {
        if ([self hasActiveVideo] && [self fullscreen] && b_videoPlayback) {
            [self hideControlsBar];
            [self.fspanel setActive];
        } else {
            [self showControlsBar];
            [self.fspanel setNonActive];
        }
    }
}

#pragma mark -
#pragma mark Fullscreen support

- (void)showFullscreenController
{
    id currentWindow = [NSApp keyWindow];
    if ([currentWindow respondsToSelector:@selector(hasActiveVideo)] && [currentWindow hasActiveVideo]) {
        if ([currentWindow respondsToSelector:@selector(fullscreen)] && [currentWindow fullscreen] && ![[currentWindow videoView] isHidden]) {

            if ([[VLCMain sharedInstance] activeVideoPlayback])
                [self.fspanel fadeIn];
        }
    }

}

@end

@implementation VLCLibraryDataSource

- (NSInteger)collectionView:(NSCollectionView *)collectionView numberOfItemsInSection:(NSInteger)section
{
    return 2;
}

- (NSCollectionViewItem *)collectionView:(NSCollectionView *)collectionView itemForRepresentedObjectAtIndexPath:(NSIndexPath *)indexPath
{
    VLCLibraryCollectionViewItem *viewItem = [collectionView makeItemWithIdentifier:VLCLibraryCellIdentifier forIndexPath:indexPath];

    viewItem.mediaTitleTextField.stringValue = @"Custom Cell Label Text";
    viewItem.mediaImageView.image = [NSImage imageNamed: @"noart.png"];

    return viewItem;
}

- (void)collectionView:(NSCollectionView *)collectionView didSelectItemsAtIndexPaths:(NSSet<NSIndexPath *> *)indexPaths
{
    NSLog(@"library selection changed: %@", indexPaths);
}

@end

@implementation VLCLibraryWindowController

- (instancetype)initWithLibraryWindow
{
    self = [super initWithWindowNibName:@"VLCLibraryWindow"];
    return self;
}

- (void)windowDidLoad
{
    VLCLibraryWindow *window = (VLCLibraryWindow *)self.window;
    [window setRestorable:NO];
    [window setExcludedFromWindowsMenu:YES];
    [window setAcceptsMouseMovedEvents:YES];
    [window setContentMinSize:NSMakeSize(f_min_window_width, f_min_window_height)];
}

@end
