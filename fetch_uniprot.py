"""Optional: refresh all nine full-JSON UniProt search results from REST.

The project already bundles the frozen 100-protein dataset, so a network
request is not needed for reproducing measured results. Running this script
may give a different cohort after UniProt updates. Requires network access.
"""

from __future__ import annotations

import json
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from build_dataset import CACHE, QUERIES


def download():
    CACHE.mkdir(parents=True, exist_ok=True)
    for family, query in QUERIES.items():
        url = 'https://rest.uniprot.org/uniprotkb/search?' + urlencode({
            'query': query, 'format': 'json', 'size': 350,
        })
        for attempt in range(5):
            try:
                req = Request(url, headers={'User-Agent': 'uniprot-guarded-codec-repro/1.0'})
                with urlopen(req, timeout=120) as response:
                    data = json.load(response)
                if 'results' not in data:
                    raise ValueError('Unexpected UniProt response')
                path = CACHE / f'{family}.json'
                path.write_text(json.dumps(data, separators=(',', ':')) + '\n')
                print(f'{family}: {len(data["results"])} entries')
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(min(2 ** attempt * 2, 16))
        time.sleep(0.5)


if __name__ == '__main__':
    download()
