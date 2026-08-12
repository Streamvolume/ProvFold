#!/usr/bin/env python3
"""Retrieve the two registered gutMDisorder exports and verify identities."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


BASE = "https://bio-computing.hrbmu.edu.cn/gutMDisorder_api/api/resource/download?fileName="
FILES = {
    "gutMDisorder_v3_Rawdata-based_Disorder_Health.xlsx": "a19bc2bd3c0b8f912fb94d4fbcc0a239a3623b345c470b1e83d2298a819d0b3c",
    "gutMDisorder_v3_Literature-based_Disorder_Health.xlsx": "8120e4b0b5b8bbdb384f374eceffb71b8ece112712a487af7a42302f2174dbb4",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    output = Path(os.environ.get("PROVFOLD_PUBLIC_INPUT_DIR", "source_inputs")).resolve()
    output.mkdir(parents=True, exist_ok=True)
    for name, expected in FILES.items():
        url = BASE + quote(name)
        request = Request(url, headers={"User-Agent": "ProvFold/0.1 source-identity-check"})
        with urlopen(request, timeout=120) as response:
            data = response.read()
        observed = digest(data)
        if observed != expected:
            raise SystemExit(
                f"HASH MISMATCH for {name}: expected {expected}, received {observed}. "
                "The official export has changed; archive it as a new source state and do not call this an exact replay."
            )
        path = output / name
        path.write_bytes(data)
        print(f"verified {name}: {observed}")


if __name__ == "__main__":
    main()
