"""Logic to handle ftp data transfer."""
from dateutil import parser
import pandas as pd
import ftplib
import io
import gzip

class FetchFromFtp:
    """This class is to retrieve the association table from the most recent GWAS Catalog release.

    Expects the ftp host address.
    - It returns the release date.
    - Fetch file or directory
    """

    def __init__(self, host: str) -> None:
        self.FTP_HOST = host

        # Initialize connection and go to folder:
        self.ftp = ftplib.FTP(self.FTP_HOST, 'anonymous', '')

    def fetch_file_list(self, path):
        # Get list of files and the date of modification:
        files = []
        self.ftp.cwd(path)
        self.ftp.dir(files.append)

        return [' '.join(x.split()[8:]) for x in files]


    def fetch_last_update_date(self, path):
        """
        This function returns the date of the most recently modified file.
        """

        # Get list of files and the date of modification:
        files = []
        self.ftp.cwd(path)
        self.ftp.dir(files.append)

        # Get all dates:
        dates = [' '.join(x.split()[5:8]) for x in files]
        dates_parsed = [parser.parse(x) for x in dates]

        release_date = max(dates_parsed)
        return release_date.strftime('%Y-%m-%d')

    def fetch_file(self, path, file):
        sio = io.BytesIO()

        def handle_binary(more_data):
            sio.write(more_data)

        self.ftp.retrbinary(f"RETR {path}/{file}", handle_binary)
        sio.seek(0)  # Go back to the start
        zippy = gzip.GzipFile(fileobj=sio)
        return zippy

    def fetch_tsv(self, path, file, skiprows=None, header='infer'):
        self.tsv_data = pd.read_csv(
            f'ftp://{self.FTP_HOST}/{path}/{file}',
            sep='\t', dtype=str, skiprows=skiprows, header=header
        )

    def close_connection(self):
        self.ftp.close()
