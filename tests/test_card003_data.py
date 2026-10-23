"""Encoding and join-identity fixtures, not CLVer observations."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import tempfile
import unittest
import zipfile
import numpy as np
import pandas as pd
from card003_data import ArchiveTables,canonical_id,numeric,source_dates,collapse_and_flag,measurement_times


class SourceHelpers(unittest.TestCase):
    def test_measurement_time_never_invents_midnight(self):
        times=measurement_times(pd.Series(["9/28/2020","9/28/2020 12:00:00 AM","9/28/2020 11:59:59 PM"]))
        self.assertTrue(pd.isna(times.iloc[0]))
        self.assertEqual(times.iloc[1],pd.Timestamp("2020-09-28 00:00:00"))
        self.assertEqual(times.iloc[2],pd.Timestamp("2020-09-28 23:59:59"))
    def test_original_encoding_and_literal_missing_codes(self):
        scratch=Path(tempfile.gettempdir()).resolve()
        with tempfile.TemporaryDirectory(dir=scratch) as directory:
            # Verify the absolute cleanup target remains inside the named test root.
            self.assertTrue(Path(directory).resolve().is_relative_to(scratch.resolve()))
            path=Path(directory)/"software_fixture.zip"
            with zipfile.ZipFile(path,"w") as z:
                z.writestr("Data Tables/Fixture.txt","PtID|Status|Value\r\n001|NA|<0.2\r\n".encode("utf-16"))
            with ArchiveTables(path) as tables:
                self.assertEqual(tables.format("Fixture.txt")["encoding"],"utf-16")
                frame=tables.read("Fixture.txt")
                self.assertEqual(frame.iloc[0].Status,"NA")
                self.assertEqual(frame.iloc[0].PtID,"001")
                self.assertEqual(frame.iloc[0].Value,"<0.2")

    def test_id_numeric_and_date_no_guessed_origin(self):
        self.assertEqual(canonical_id("001.0"),"1")
        with self.assertRaises(ValueError): canonical_id("1.5")
        with self.assertRaises(ValueError): source_dates(["12345"])
        self.assertEqual(source_dates(["2024-01-03"]).iloc[0],pd.Timestamp("2024-01-03"))
        parsed=numeric(pd.Series(["0.2","<0.2","inf",""]))
        self.assertEqual(parsed.iloc[0],.2)
        self.assertTrue(parsed.iloc[1:].isna().all())

    def test_conflicts_exclude_all_rows_in_key(self):
        f=pd.DataFrame(dict(PtID=[1,1,2,2,3],Visit=[13]*5,Value=[.2,.2,.3,.4,.1],RecID=[1,2,3,4,5]))
        kept,bad,audit=collapse_and_flag(f,["PtID","Visit"],["Value"])
        self.assertEqual(set(kept.PtID),{1,3})
        self.assertEqual(set(bad.PtID),{2})
        self.assertEqual(audit["analytical_duplicate_rows"],1)
        self.assertEqual(audit["conflicting_keys"],1)

if __name__=="__main__": unittest.main()
