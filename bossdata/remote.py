# -*- coding: utf-8 -*-

"""Download BOSS data files from a remote server.
"""

from __future__ import division,print_function

import os
import os.path

import requests

class Manager(object):
    """Manage downloads of BOSS data via HTTP.

    The default mapping from remote to local filenames is to mirror the remote file hierarchy
    on the local disk.  The normal mode of operation is to establish the local root for the
    mirror using the BOSSDATA_ROOT environment variable.

    Args:
        release(str): Tag describing the release to download from, which determines the
            base URL of all BOSS data files. Currently only 'DR12' is supported.
        local_root(str): Local path to use as the root of the locally mirrored file
            hierarchy. If this arg is None, then the value of the BOSSDATA_ROOT environment
            variable, if any, will be used instead.  If a value is provided, it should identify
            an existing writeable directory.

    Raises:
        RuntimeError: No such directory local_root.
    """
    def __init__(self,release = 'DR12',local_root = None):

        if release == 'DR12':
            self.url_prefix = 'http://dr12.sdss3.org'
        else:
            raise ValueError('Invalid release {}.'.format(source))

        self.local_root = local_root
        if self.local_root is None:
            self.local_root = os.getenv('BOSSDATA_ROOT')
        if self.local_root is not None and not os.path.isdir(self.local_root):
            raise RuntimeError('Cannot use non-existent path {} as local root.'.format(self.local_root))

    def download(self,remote_path,local_path = None,chunk_size = 4096):
        """Download a single BOSS data file.

        Downloads are streamed so that the memory requirements are independent of the file size.

        TODO: automatic decompression, optional progress feedback, check for enough disk space.

        Args:
            remote_path(str): The full path to the remote file relative to the remote server root.
            local_path(str): The (absolute or relative) path of the local file to write. When no
                local_path is provided, the default behavior is to mirror the remote file hierarchy
                under the local_root, creating new directories as needed.
            chunk_size(int): Size of data chunks to use for the streaming download. Larger sizes
                will potentially download faster but also require more memory.

        Returns:
            str: Absolute local path of the downloaded file.

        Raises:
            RuntimeError: no local_path provided and no local_root specified when this Manager
                was created.
        """
        # Prepare the HTTP request.
        request = requests.get(self.url_prefix + remote_path,stream = True)

        # Mirror the remote directory layout under the local root by default.
        if local_path is None and self.local_root is not None:
            local_path = os.path.join(self.local_root,remote_path)

        if local_path is None:
            raise RuntimeError('Cannot mirror without a local root (try setting BOSSDATA_ROOT).')

        # Create local directories as needed.
        local_path = os.path.abspath(local_path)
        parent_path = os.path.dirname(local_path)
        if not os.path.exists(parent_path):
            os.makedirs(parent_path)

        # Stream the request response binary content into the local file.
        with open(local_path,'wb') as f:
            for chunk in request.iter_content(chunk_size = chunk_size):
                f.write(chunk)

        return local_path
