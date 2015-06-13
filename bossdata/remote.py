# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Download BOSS data files from a remote server.
"""

from __future__ import division, print_function

import os
import stat
import os.path
import math

import requests

from progressbar import ProgressBar, Percentage, Bar, FileTransferSpeed


class Manager(object):
    """Manage downloads of BOSS data via HTTP.

    The default mapping from remote to local filenames is to mirror the remote file hierarchy
    on the local disk.  The normal mode of operation is to establish the local root for the
    mirror using the BOSS_LOCAL_ROOT environment variable. When the constructor is called
    with no arguments, it will raise a ValueError if either BOSS_DATA_URL or BOSS_LOCAL_ROOT
    is not set.

    Args:
        data_url(str): Base URL of all BOSS data files. A trailing / on the URL is optional. If
            this arg is None, then the value of the BOSS_DATA_URL environment variable we be
            used instead.
        local_root(str): Local path to use as the root of the locally mirrored file
            hierarchy. If this arg is None, then the value of the BOSS_LOCAL_ROOT environment
            variable, if any, will be used instead.  If a value is provided, it should identify
            an existing writeable directory.

    Raises:
        ValueError: No such directory local_root or missing data_url.
    """
    def __init__(self, data_url=None, local_root=None):

        self.data_url = data_url
        if self.data_url is None:
            self.data_url = os.getenv('BOSS_DATA_URL')
        if self.data_url is None:
            raise ValueError('No data URL specified: try setting $BOSS_DATA_URL.')
        self.data_url = self.data_url.rstrip('/')

        self.local_root = local_root
        if self.local_root is None:
            self.local_root = os.getenv('BOSS_LOCAL_ROOT')
            if self.local_root is None:
                raise ValueError('No local root specified: try setting $BOSS_LOCAL_ROOT.')
        if not os.path.isdir(self.local_root):
            raise ValueError('Cannot use non-existent path "{}" as local root.'.format(
                self.local_root))

    def download(self, remote_path, local_path, chunk_size=4096, progress_min_size=10):
        """Download a single BOSS data file.

        Downloads are streamed so that the memory requirements are independent of the
        file size. During the download, the file is written to its final location but
        with '.downloading' appended to the file name. This means than any download
        that is interrupted or fails will normally not lead to an incomplete file
        being returned by a subsequent call to :meth:`get`. Instead, the file will
        be re-downloaded. Tere is no facility for resuming a previous partial download.
        After a successful download, the file is renamed to its final location and
        has its permission bits set to read only (to prevent accidental modifications
        of files that are supposed to exactly mirror the remote file system).

        Args:
            remote_path(str): The full path to the remote file relative to the remote
                server root, which should normally be obtained using :class:`bossdata.path`
                methods.
            local_path(str): The (absolute or relative) path of the local file to write.
            chunk_size(int): Size of data chunks to use for the streaming download. Larger
                sizes will potentially download faster but also require more memory.
            progress_min_size(int): Display a text progress bar for any downloads whose size
                in Mb exceeds this value.

        Returns:
            str: Absolute local path of the downloaded file.

        Raises:
            ValueError: local_path directory does not exist.
            RuntimeError: HTTP request returned an error status.
        """
        if not local_path:
            raise ValueError('Missing required argument local_path.')
        local_path = os.path.abspath(local_path)

        # Check that the local path points to an existing directory.
        if not os.path.isdir(os.path.dirname(local_path)):
            raise ValueError('local_path directory does not exist: {}.'.format(
                os.path.dirname(local_path)))

        # Prepare the HTTP request. For details on the timeout parameter see
        # http://docs.python-requests.org/en/latest/user/advanced/#timeouts
        url = self.data_url + '/' + remote_path.lstrip('/')
        try:
            request = requests.get(url, stream=True, timeout=(3.05, 27))
            if request.status_code != requests.codes.ok:
                raise RuntimeError('HTTP request returned error code {0}.'.format(
                    request.status_code))
        except requests.exceptions.RequestException as e:
            raise RuntimeError('HTTP request failed: {}.'.format(str(e)))

        # Check that there is enough free space, if possible.
        progress_bar = None
        file_size = request.headers.get('content-length', None)
        if file_size is not None:
            file_size = int(file_size)
            parent_path = os.path.dirname(local_path)
            file_stat = os.statvfs(parent_path)
            free_space = file_stat.f_bavail * file_stat.f_frsize
            Mb = 1 << 20
            if file_size + 1 * Mb > free_space:
                raise RuntimeError('File size ({:.1f}Mb) exceeds free space for {}.'.format(
                    file_size / (1.0 * Mb), local_path))
            if file_size > progress_min_size * Mb:
                label = os.path.basename(local_path)
                progress_bar = ProgressBar(
                    widgets=[label, ' ', Percentage(), Bar(), ' ', FileTransferSpeed()],
                    maxval=math.ceil(file_size / chunk_size)).start()

        # Stream the request response binary content into a temporary file.
        progress = 0
        with open(local_path + '.downloading', 'wb') as f:
            for chunk in request.iter_content(chunk_size=chunk_size):
                f.write(chunk)
                if progress_bar:
                    progress += 1
                    progress_bar.update(progress)
        # Make the temporary file read only by anyone.
        os.chmod(local_path + '.downloading', stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)
        # Move the temporary file to its permanent location.
        os.rename(local_path + '.downloading', local_path)
        if progress_bar:
            progress_bar.finish()

        return local_path

    def local_path(self, remote_path):
        """Get the local path corresponding to a remote path.

        Does not check that the file or its parent directory exists. Use :meth:`get` to
        ensure that the file exists, downloading it if necessary.

        Args:
            remote_path(str): The full path to the remote file relative to the remote
                server root, which should normally be obtained using :class:`bossdata.path`
                methods.

        Returns:
            str: Absolute local path of the local file that mirrors the remote file.

        Raises:
            RuntimeError: No local_root specified when this manager was created.
        """
        if self.local_root is None:
            raise RuntimeError('No local root specified (try setting BOSS_LOCAL_ROOT).')
        return os.path.abspath(os.path.join(self.local_root, remote_path.lstrip('/')))

    def get(self, remote_path):
        """Get a local file that mirrors a remote file, downloading the file if necessary.

        Args:
            remote_path(str): The full path to the remote file relative to the remote
                server root, which should normally be obtained using :class:`bossdata.path`
                methods.

        Returns:
            str: Absolute local path of the local file that mirrors the remote file.
        """
        local_path = self.local_path(remote_path)
        if os.path.isfile(local_path):
            return local_path

        # If we get here, the file is not available locally so try to download it now.
        # Create local directories as needed.
        parent_path = os.path.dirname(local_path)
        if not os.path.exists(parent_path):
            # There is a potential race condition if other processes are running.
            # In python >= 3.2 we would use the new exist_ok=True option, but here we instead
            # silently ignore a "FileExists" OSError (errno 17).
            try:
                os.makedirs(parent_path)
            except OSError as e:
                if e.errno != 17:
                    raise e
        return self.download(remote_path, local_path)
