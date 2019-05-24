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
#import "extensions/NSFont+VLCAdditions.h"
#import "extensions/NSColor+VLCAdditions.h"
#import "extensions/NSView+VLCAdditions.h"
#import "main/VLCMain.h"

#import "playlist/VLCPlaylistController.h"
#import "playlist/VLCPlaylistDataSource.h"

#import "library/VLCLibraryController.h"
#import "library/VLCLibraryAudioDataSource.h"
#import "library/VLCLibraryVideoDataSource.h"
#import "library/VLCLibraryCollectionViewItem.h"
#import "library/VLCLibraryModel.h"
#import "library/VLCLibraryCollectionViewSupplementaryElementView.h"
#import "library/VLCLibraryAlternativeAudioViewController.h"

#import "media-source/VLCMediaSourceCollectionViewItem.h"
#import "media-source/VLCMediaSourceDataSource.h"

#import "windows/mainwindow/VLCControlsBarCommon.h"
#import "windows/video/VLCFSPanelController.h"
#import "windows/video/VLCVoutView.h"
#import "windows/video/VLCVideoOutputProvider.h"

const CGFloat VLCLibraryWindowMinimalWidth = 604.;
const CGFloat VLCLibraryWindowMinimalHeight = 307.;
const CGFloat VLCLibraryWindowPlaylistRowHeight = 72.;
const CGFloat VLCLibraryWindowSmallRowHeight = 24.;
const CGFloat VLCLibraryWindowLargeRowHeight = 50.;

@interface VLCLibraryWindow ()
{
    VLCPlaylistDataSource *_playlistDataSource;
    VLCLibraryVideoDataSource *_libraryVideoDataSource;
    VLCLibraryAudioDataSource *_libraryAudioDataSource;
    VLCLibraryGroupDataSource *_libraryAudioGroupDataSource;
    VLCMediaSourceDataSource *_mediaSourceDataSource;
    VLCLibraryAlternativeAudioViewController *_alternativeAudioViewController;

    VLCPlaylistController *_playlistController;

    NSRect _windowFrameBeforePlayback;

    VLCFSPanelController *_fspanel;
}
@end

@implementation VLCLibraryWindow

- (void)awakeFromNib
{
    VLCMain *mainInstance = [VLCMain sharedInstance];
    _playlistController = [mainInstance playlistController];

    self.videoView = [[VLCVoutView alloc] initWithFrame:self.mainSplitView.frame];
    self.videoView.hidden = YES;
    self.videoView.translatesAutoresizingMaskIntoConstraints = NO;
    [self.contentView addSubview:self.videoView];
    [self.contentView addConstraint:[NSLayoutConstraint constraintWithItem:self.videoView attribute:NSLayoutAttributeWidth relatedBy:NSLayoutRelationEqual toItem:self.mainSplitView attribute:NSLayoutAttributeWidth multiplier:1. constant:1.]];
    [self.contentView addConstraint:[NSLayoutConstraint constraintWithItem:self.videoView attribute:NSLayoutAttributeHeight relatedBy:NSLayoutRelationEqual toItem:self.mainSplitView attribute:NSLayoutAttributeHeight multiplier:1. constant:1.]];
    [self.contentView addConstraint:[NSLayoutConstraint constraintWithItem:self.videoView attribute:NSLayoutAttributeCenterX relatedBy:NSLayoutRelationEqual toItem:self.mainSplitView attribute:NSLayoutAttributeCenterX multiplier:1. constant:1.]];
    [self.contentView addConstraint:[NSLayoutConstraint constraintWithItem:self.videoView attribute:NSLayoutAttributeCenterY relatedBy:NSLayoutRelationEqual toItem:self.mainSplitView attribute:NSLayoutAttributeCenterY multiplier:1. constant:1.]];

    NSNotificationCenter *notificationCenter = [NSNotificationCenter defaultCenter];
    [notificationCenter addObserver:self
                           selector:@selector(shouldShowFullscreenController:)
                               name:VLCVideoWindowShouldShowFullscreenController
                             object:nil];
    [notificationCenter addObserver:self
                           selector:@selector(updateLibraryRepresentation:)
                               name:VLCLibraryModelAudioMediaListUpdated
                             object:nil];
    [notificationCenter addObserver:self
                           selector:@selector(updateLibraryRepresentation:)
                               name:VLCLibraryModelVideoMediaListUpdated
                             object:nil];
    [notificationCenter addObserver:self
                           selector:@selector(updateLibraryRepresentation:)
                               name:VLCLibraryModelRecentMediaListUpdated
                             object:nil];
    [notificationCenter addObserver:self
                           selector:@selector(shuffleStateUpdated:)
                               name:VLCPlaybackOrderChanged
                             object:nil];
    [notificationCenter addObserver:self
                           selector:@selector(repeatStateUpdated:)
                               name:VLCPlaybackRepeatChanged
                             object:nil];

    if (@available(macOS 10_14, *)) {
        [[NSApplication sharedApplication] addObserver:self
                                            forKeyPath:@"effectiveAppearance"
                                               options:0
                                               context:nil];
    }

    _fspanel = [[VLCFSPanelController alloc] init];
    [_fspanel showWindow:self];

    _segmentedTitleControl.segmentCount = 5;
    [_segmentedTitleControl setTarget:self];
    [_segmentedTitleControl setAction:@selector(segmentedControlAction:)];
    [_segmentedTitleControl setLabel:_NS("Video") forSegment:0];
    [_segmentedTitleControl setLabel:_NS("Music") forSegment:1];
    [_segmentedTitleControl setLabel:_NS("Music") forSegment:2];
    [_segmentedTitleControl setLabel:_NS("Local Network") forSegment:3];
    [_segmentedTitleControl setLabel:_NS("Internet") forSegment:4];
    [_segmentedTitleControl sizeToFit];
    [_segmentedTitleControl setSelectedSegment:0];

    _playlistDataSource = [[VLCPlaylistDataSource alloc] init];
    _playlistDataSource.playlistController = _playlistController;
    _playlistDataSource.tableView = _playlistTableView;
    _playlistController.playlistDataSource = _playlistDataSource;

    _playlistTableView.dataSource = _playlistDataSource;
    _playlistTableView.delegate = _playlistDataSource;
    _playlistTableView.rowHeight = VLCLibraryWindowPlaylistRowHeight;
    [_playlistTableView reloadData];

    _libraryVideoDataSource = [[VLCLibraryVideoDataSource alloc] init];
    _libraryVideoDataSource.libraryModel = mainInstance.libraryController.libraryModel;
    _libraryVideoDataSource.recentMediaCollectionView = _recentVideoLibraryCollectionView;
    _libraryVideoDataSource.libraryMediaCollectionView = _videoLibraryCollectionView;
    _videoLibraryCollectionView.dataSource = _libraryVideoDataSource;
    _videoLibraryCollectionView.delegate = _libraryVideoDataSource;
    [_videoLibraryCollectionView registerClass:[VLCLibraryCollectionViewItem class] forItemWithIdentifier:VLCLibraryCellIdentifier];
    [_videoLibraryCollectionView registerClass:[VLCLibraryCollectionViewSupplementaryElementView class]
               forSupplementaryViewOfKind:NSCollectionElementKindSectionHeader
                           withIdentifier:VLCLibrarySupplementaryElementViewIdentifier];
    [(NSCollectionViewFlowLayout *)_videoLibraryCollectionView.collectionViewLayout setHeaderReferenceSize:[VLCLibraryCollectionViewSupplementaryElementView defaultHeaderSize]];
    _recentVideoLibraryCollectionView.dataSource = _libraryVideoDataSource;
    _recentVideoLibraryCollectionView.delegate = _libraryVideoDataSource;
    [_recentVideoLibraryCollectionView registerClass:[VLCLibraryCollectionViewItem class] forItemWithIdentifier:VLCLibraryCellIdentifier];

    _libraryAudioDataSource = [[VLCLibraryAudioDataSource alloc] init];
    _libraryAudioDataSource.libraryModel = mainInstance.libraryController.libraryModel;
    _libraryAudioDataSource.categorySelectionTableView = _audioCategorySelectionTableView;
    _libraryAudioDataSource.collectionSelectionTableView = _audioCollectionSelectionTableView;
    _libraryAudioDataSource.groupSelectionTableView = _audioGroupSelectionTableView;
    _audioCategorySelectionTableView.dataSource = _libraryAudioDataSource;
    _audioCategorySelectionTableView.delegate = _libraryAudioDataSource;
    _audioCategorySelectionTableView.rowHeight = VLCLibraryWindowSmallRowHeight;
    _audioCollectionSelectionTableView.dataSource = _libraryAudioDataSource;
    _audioCollectionSelectionTableView.delegate = _libraryAudioDataSource;
    _audioCollectionSelectionTableView.rowHeight = VLCLibraryWindowLargeRowHeight;
    _libraryAudioGroupDataSource = [[VLCLibraryGroupDataSource alloc] init];
    _libraryAudioDataSource.groupDataSource = _libraryAudioGroupDataSource;
    _audioGroupSelectionTableView.dataSource = _libraryAudioGroupDataSource;
    _audioGroupSelectionTableView.delegate = _libraryAudioGroupDataSource;
    _audioGroupSelectionTableView.rowHeight = 450.;

    _mediaSourceDataSource = [[VLCMediaSourceDataSource alloc] init];
    _mediaSourceDataSource.collectionView = _mediaSourceCollectionView;
    _mediaSourceCollectionView.dataSource = _mediaSourceDataSource;
    _mediaSourceCollectionView.delegate = _mediaSourceDataSource;
    [_mediaSourceCollectionView registerClass:[VLCMediaSourceCollectionViewItem class] forItemWithIdentifier:VLCMediaSourceCellIdentifier];

    self.upNextLabel.font = [NSFont VLClibrarySectionHeaderFont];
    self.upNextLabel.stringValue = _NS("Up next");
    NSAttributedString *attributedTitle = [[NSAttributedString alloc] initWithString:_NS("Clear queue")
                                                                          attributes:@{NSFontAttributeName : [NSFont VLClibraryButtonFont],
                                                                                       NSForegroundColorAttributeName : [NSColor VLClibraryHighlightColor]}];
    self.clearPlaylistButton.attributedTitle = attributedTitle;
    [self updateColorsBasedOnAppearance];

    _alternativeAudioViewController = [[VLCLibraryAlternativeAudioViewController alloc] init];
    _alternativeAudioViewController.collectionView = self.alternativeAudioCollectionView;
    _alternativeAudioViewController.segmentedControl = self.alternativeAudioSegmentedControl;
    _alternativeAudioViewController.libraryModel = mainInstance.libraryController.libraryModel;
    [_alternativeAudioViewController setupAppearance];

    [self segmentedControlAction:nil];
    [self repeatStateUpdated:nil];
    [self shuffleStateUpdated:nil];
}

- (void)dealloc
{
    [[NSNotificationCenter defaultCenter] removeObserver:self];
    if (@available(macOS 10_14, *)) {
        [[NSApplication sharedApplication] removeObserver:self forKeyPath:@"effectiveAppearance"];
    }
}

- (void)observeValueForKeyPath:(NSString *)keyPath
                      ofObject:(id)object
                        change:(NSDictionary<NSKeyValueChangeKey,id> *)change
                       context:(void *)context
{
    [self updateColorsBasedOnAppearance];
}

- (void)updateColorsBasedOnAppearance
{
    if (self.contentView.shouldShowDarkAppearance) {
        self.upNextLabel.textColor = [NSColor VLClibraryDarkTitleColor];
        self.upNextSeparator.borderColor = [NSColor VLClibrarySeparatorDarkColor];
        self.clearPlaylistSeparator.borderColor = [NSColor VLClibrarySeparatorDarkColor];
    } else {
        self.upNextLabel.textColor = [NSColor VLClibraryLightTitleColor];
        self.upNextSeparator.borderColor = [NSColor VLClibrarySeparatorLightColor];
        self.clearPlaylistSeparator.borderColor = [NSColor VLClibrarySeparatorLightColor];
    }
}

#pragma mark - playmode state display and interaction

- (IBAction)shuffleAction:(id)sender
{
    if (_playlistController.playbackOrder == VLC_PLAYLIST_PLAYBACK_ORDER_NORMAL) {
        _playlistController.playbackOrder = VLC_PLAYLIST_PLAYBACK_ORDER_RANDOM;
    } else {
        _playlistController.playbackOrder = VLC_PLAYLIST_PLAYBACK_ORDER_NORMAL;
    }
}

- (void)shuffleStateUpdated:(NSNotification *)aNotification
{
    if (_playlistController.playbackOrder == VLC_PLAYLIST_PLAYBACK_ORDER_NORMAL) {
        self.shufflePlaylistButton.image = [NSImage imageNamed:@"shuffleOff"];
    } else {
        self.shufflePlaylistButton.image = [NSImage imageNamed:@"shuffleOn"];
    }
}

- (IBAction)repeatAction:(id)sender
{
    enum vlc_playlist_playback_repeat currentRepeatState = _playlistController.playbackRepeat;
    switch (currentRepeatState) {
        case VLC_PLAYLIST_PLAYBACK_REPEAT_ALL:
            _playlistController.playbackRepeat = VLC_PLAYLIST_PLAYBACK_REPEAT_NONE;
            break;
        case VLC_PLAYLIST_PLAYBACK_REPEAT_CURRENT:
            _playlistController.playbackRepeat = VLC_PLAYLIST_PLAYBACK_REPEAT_ALL;
            break;

        default:
            _playlistController.playbackRepeat = VLC_PLAYLIST_PLAYBACK_REPEAT_CURRENT;
            break;
    }
}

- (void)repeatStateUpdated:(NSNotification *)aNotification
{
    enum vlc_playlist_playback_repeat currentRepeatState = _playlistController.playbackRepeat;
    switch (currentRepeatState) {
        case VLC_PLAYLIST_PLAYBACK_REPEAT_ALL:
            self.repeatPlaylistButton.image = [NSImage imageNamed:@"repeatAll"];
            break;
        case VLC_PLAYLIST_PLAYBACK_REPEAT_CURRENT:
            self.repeatPlaylistButton.image = [NSImage imageNamed:@"repeatOne"];
            break;

        default:
            self.repeatPlaylistButton.image = [NSImage imageNamed:@"repeatOff"];
            break;
    }
}

#pragma mark - misc. user interactions

- (void)segmentedControlAction:(id)sender
{
    switch (_segmentedTitleControl.selectedSegment) {
        case 0:
            if (_mediaSourceScrollView.superview != nil) {
                [_mediaSourceScrollView removeFromSuperview];
            }
            if (_audioLibrarySplitView.superview != nil) {
                [_audioLibrarySplitView removeFromSuperview];
            }
            if (_alternativeAudioView.superview != nil) {
                [_alternativeAudioView removeFromSuperview];
            }
            if (_videoLibraryStackView.superview == nil) {
                _videoLibraryStackView.translatesAutoresizingMaskIntoConstraints = NO;
                [_libraryTargetView addSubview:_videoLibraryStackView];
                NSDictionary *dict = NSDictionaryOfVariableBindings(_videoLibraryStackView);
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[_videoLibraryStackView(>=572.)]|" options:0 metrics:0 views:dict]];
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"V:|[_videoLibraryStackView(>=444.)]|" options:0 metrics:0 views:dict]];
            }
            [_videoLibraryCollectionView reloadData];
            [_recentVideoLibraryCollectionView reloadData];
            break;

        case 1:
            if (_mediaSourceScrollView.superview != nil) {
                [_mediaSourceScrollView removeFromSuperview];
            }
            if (_videoLibraryStackView.superview != nil) {
                [_videoLibraryStackView removeFromSuperview];
            }
            if (_alternativeAudioView.superview != nil) {
                [_alternativeAudioView removeFromSuperview];
            }
            if (_audioLibrarySplitView.superview == nil) {
                _audioLibrarySplitView.translatesAutoresizingMaskIntoConstraints = NO;
                [_libraryTargetView addSubview:_audioLibrarySplitView];
                NSDictionary *dict = NSDictionaryOfVariableBindings(_audioLibrarySplitView);
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[_audioLibrarySplitView(>=572.)]|" options:0 metrics:0 views:dict]];
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"V:|[_audioLibrarySplitView(>=444.)]|" options:0 metrics:0 views:dict]];
            }
            [_audioCategorySelectionTableView reloadData];
            [_audioCollectionSelectionTableView reloadData];
            break;

        case 2:
            if (_mediaSourceScrollView.superview != nil) {
                [_mediaSourceScrollView removeFromSuperview];
            }
            if (_videoLibraryStackView.superview != nil) {
                [_videoLibraryStackView removeFromSuperview];
            }
            if (_audioLibrarySplitView.superview != nil) {
                [_audioLibrarySplitView removeFromSuperview];
            }
            if (_alternativeAudioView.superview == nil) {
                _alternativeAudioView.translatesAutoresizingMaskIntoConstraints = NO;
                [_libraryTargetView addSubview:_alternativeAudioView];
                NSDictionary *dict = NSDictionaryOfVariableBindings(_alternativeAudioView);
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[_alternativeAudioView(>=572.)]|" options:0 metrics:0 views:dict]];
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"V:|[_alternativeAudioView(>=444.)]|" options:0 metrics:0 views:dict]];
            }
            [_alternativeAudioViewController reloadAppearance];
            break;

        default:
            if (_videoLibraryStackView.superview != nil) {
                [_videoLibraryStackView removeFromSuperview];
            }
            if (_audioLibrarySplitView.superview != nil) {
                [_audioLibrarySplitView removeFromSuperview];
            }
            if (_alternativeAudioView.superview != nil) {
                [_alternativeAudioView removeFromSuperview];
            }
            if (_mediaSourceScrollView.superview == nil) {
                _mediaSourceScrollView.translatesAutoresizingMaskIntoConstraints = NO;
                [_libraryTargetView addSubview:_mediaSourceScrollView];
                NSDictionary *dict = NSDictionaryOfVariableBindings(_mediaSourceScrollView);
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"H:|[_mediaSourceScrollView(>=572.)]|" options:0 metrics:0 views:dict]];
                [_libraryTargetView addConstraints:[NSLayoutConstraint constraintsWithVisualFormat:@"V:|[_mediaSourceScrollView(>=444.)]|" options:0 metrics:0 views:dict]];
            }
            [_mediaSourceDataSource loadMediaSources];
            [_mediaSourceCollectionView reloadData];
            break;
    }
}

- (IBAction)playlistDoubleClickAction:(id)sender
{
    NSInteger selectedRow = self.playlistTableView.selectedRow;
    if (selectedRow == -1)
        return;

    [[[VLCMain sharedInstance] playlistController] playItemAtIndex:selectedRow];
}

- (IBAction)clearPlaylist:(id)sender
{
    [_playlistController clearPlaylist];
}

#pragma mark - video output controlling

- (void)videoPlaybackWillBeStarted
{
    if (!self.fullscreen)
        _windowFrameBeforePlayback = [self frame];
}

- (void)enableVideoPlaybackAppearance
{
    [self.videoView setHidden:NO];

    if (self.nativeFullscreenMode) {
        if ([self hasActiveVideo] && [self fullscreen]) {
            [self hideControlsBar];
            [_fspanel shouldBecomeActive:nil];
        }
    }
}

- (void)disableVideoPlaybackAppearance
{
    if (!self.nonembedded
        && (!self.nativeFullscreenMode || (self.nativeFullscreenMode && !self.fullscreen))
        && _windowFrameBeforePlayback.size.width > 0
        && _windowFrameBeforePlayback.size.height > 0) {

        // only resize back to minimum view of this is still desired final state
        CGFloat f_threshold_height = VLCVideoWindowCommonMinimalHeight + [self.controlsBar height];
        if (_windowFrameBeforePlayback.size.height > f_threshold_height) {
            if ([[VLCMain sharedInstance] isTerminating]) {
                [self setFrame:_windowFrameBeforePlayback display:YES];
            } else {
                [[self animator] setFrame:_windowFrameBeforePlayback display:YES];
            }
        }
    }

    _windowFrameBeforePlayback = NSMakeRect(0, 0, 0, 0);

    [self makeFirstResponder: _playlistTableView];
    [[[VLCMain sharedInstance] voutProvider] updateWindowLevelForHelperWindows: NSNormalWindowLevel];

    // restore alpha value to 1 for the case that macosx-opaqueness is set to < 1
    [self setAlphaValue:1.0];
    [self.videoView setHidden:YES];

    if (self.nativeFullscreenMode) {
        [self showControlsBar];
        [_fspanel shouldBecomeInactive:nil];
    }
}

#pragma mark - library representation and interaction
- (void)updateLibraryRepresentation:(NSNotification *)aNotification
{
    [_videoLibraryCollectionView reloadData];
    [_recentVideoLibraryCollectionView reloadData];
}

#pragma mark -
#pragma mark Fullscreen support

- (void)shouldShowFullscreenController:(NSNotification *)aNotification
{
    id currentWindow = [NSApp keyWindow];
    if ([currentWindow respondsToSelector:@selector(hasActiveVideo)] && [currentWindow hasActiveVideo]) {
        if ([currentWindow respondsToSelector:@selector(fullscreen)] && [currentWindow fullscreen] && ![[currentWindow videoView] isHidden]) {
            if ([[VLCMain sharedInstance] activeVideoPlayback]) {
                [_fspanel fadeIn];
            }
        }
    }

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
    [window setContentMinSize:NSMakeSize(VLCLibraryWindowMinimalWidth, VLCLibraryWindowMinimalHeight)];
}

@end
