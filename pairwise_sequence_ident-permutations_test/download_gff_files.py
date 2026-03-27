#!/usr/bin/env python3
"""
Download GFF files from Google Drive folder to /tmp/gff_files
"""

import os
import gdown

GDRIVE_FOLDER_URL = "https://drive.google.com/drive/folders/19z0_k9jxSYfqiD112sbSOcp6PiOAlkti"
GDRIVE_FOLDER_ID  = "19z0_k9jxSYfqiD112sbSOcp6PiOAlkti"
OUTPUT_DIR = "gff_files"

def main():
    print("=" * 60)
    print("Downloading GFF files from Google Drive")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Download entire folder
    print(f"\nDownloading folder: {GDRIVE_FOLDER_ID}")
    print(f"Destination: {OUTPUT_DIR}")

    gdown.download_folder(
        id=GDRIVE_FOLDER_ID,
        output=OUTPUT_DIR,
        quiet=False,
        use_cookies=False
    )

    # List downloaded files
    files = os.listdir(OUTPUT_DIR)
    print(f"\n{'=' * 60}")
    print(f"Downloaded {len(files)} files:")
    for f in sorted(files):
        print(f"  - {f}")
    print(f"{'=' * 60}")

if __name__ == "__main__":
    main()
