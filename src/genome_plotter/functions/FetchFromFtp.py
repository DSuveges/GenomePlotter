"""Module for fetching data from FTP servers."""

from __future__ import annotations

import ftplib
import gzip
import io
from typing import Any

import pandas as pd
from dateutil import parser


class FetchFromFtp:
    """Class to retrieve data from FTP servers.

    This class retrieves the association table from the most recent GWAS
    Catalog release. It expects the FTP host address and returns the
    release date.
    """

    def __init__(self: FetchFromFtp, url: str) -> None:
        """Initialize FTP connection.

        Args:
            url (str): FTP host address.
        """
        self.FTP_HOST = url

        # Initialize connection and go to folder:
        self.ftp = ftplib.FTP(self.FTP_HOST, "anonymous", "")

    def fetch_file_list(self: FetchFromFtp, path: str) -> list[str]:
        """Extract the list of files in the specified path.

        Args:
            path (str): The path to the directory.

        Returns:
            list[str]: The list of files in the directory.

        Raises:
            ValueError: If no files are found in the specified path.
        """
        # Get list of files and the date of modification:
        files: list[str] = []
        self.ftp.cwd(path)
        self.ftp.dir(files.append)

        files = [" ".join(x.split()[8:]) for x in files]

        if files is None or len(files) == 0:
            raise ValueError("No files found in the specified path.")

        return files

    def fetch_last_update_date(self: FetchFromFtp, path: str) -> str:
        """Return the date of the most recently modified file.

        Args:
            path (str): Directory path on FTP server.

        Returns:
            str: Date string in YYYY-MM-DD format.
        """
        # Get list of files and the date of modification:
        files: list[str] = []
        self.ftp.cwd(path)
        self.ftp.dir(files.append)

        # Get all dates:
        dates = [" ".join(x.split()[5:8]) for x in files]
        dates_parsed = [parser.parse(x) for x in dates]

        release_date = max(dates_parsed)
        return release_date.strftime("%Y-%m-%d")

    def fetch_file(self: FetchFromFtp, path: str, file: str) -> gzip.GzipFile:
        """Fetch and decompress a gzipped file from FTP.

        Args:
            path (str): Directory path on FTP server.
            file (str): File name to fetch.

        Returns:
            gzip.GzipFile: Decompressed file object.
        """
        sio = io.BytesIO()

        def handle_binary(more_data: bytes) -> None:
            sio.write(more_data)

        self.ftp.retrbinary(f"RETR {path}/{file}", handle_binary)
        sio.seek(0)  # Go back to the start
        zippy = gzip.GzipFile(fileobj=sio)
        return zippy

    def fetch_tsv(
        self: FetchFromFtp,
        path: str,
        file: str,
        skiprows: int | None = None,
        header: Any = "infer",
    ) -> None:
        """Fetch a TSV file from FTP and store as DataFrame.

        Args:
            path (str): Directory path on FTP server.
            file (str): File name to fetch.
            skiprows (int | None): Number of rows to skip.
            header (Any): Row to use as header.
        """
        self.tsv_data = pd.read_csv(
            f"ftp://{self.FTP_HOST}/{path}/{file}",
            sep="\t",
            dtype=str,
            skiprows=skiprows,
            header=header,
        )

    def close_connection(self: FetchFromFtp) -> None:
        """Close the FTP connection."""
        self.ftp.close()
