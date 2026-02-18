#!/usr/bin/env python3

import pandas as pd

if __name__ == "__main__":
    print ("This is a class file, import me into a script")

class FileImporter:

    """A class to import qPCR data files formatted specifically from the QuantStudio 7.

    Usage:
      importer = FileImporter("path/to/file.xlsx")
      df = importer.import_file()
    """

    def __init__(self, filename: str):
        if not filename:
            raise ValueError("filename must be provided to FileImporter")
        self.filename = filename

    def import_file(self) -> pd.DataFrame:
        """Import the excel file and return a pandas DataFrame.

        The method will read the Results sheet and return a tidy DataFrame with
        Sample Name, Target Name, and CT. It will normalize common 'undetermined'
        values to NaN and drop rows with missing CT.
        """
        try:
            df = pd.read_excel(
                self.filename,
                sheet_name="Results",
                skiprows=46,
                usecols=["Sample Name", "Target Name", "CT", "Task"],
                na_values=["Undetermined", "NTC"]
            )
        except Exception as e:
            raise RuntimeError(f"failed to read Excel file '{self.filename}': {e}")

        # Keep rows with CT present and drop fully empty rows
        df = df.dropna(subset=["Sample Name", "CT"])
        # Reset index for cleanliness
        df = df.reset_index(drop=True)
        return df
    
