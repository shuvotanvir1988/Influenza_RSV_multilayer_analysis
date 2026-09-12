#!/usr/bin/env python3

from pathlib import Path
from urllib.request import urlretrieve

ROOT = Path.home() / "Influenza_RSV_Project"

OUT = ROOT / "data" / "raw" / "DS001_GSE38900" / "platform"

OUT.mkdir(parents=True, exist_ok=True)

urls = {
    "GPL10558":
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL10558&form=text&view=full",

    "GPL6884":
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL6884&form=text&view=full"
}

for gpl, url in urls.items():

    outfile = OUT / f"{gpl}.txt"

    print(f"Downloading {gpl}")

    urlretrieve(url, outfile)

    print(outfile)

print("\nDone.")
