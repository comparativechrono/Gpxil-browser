"""Static file server with HTTP range requests, for trying the browser locally.

    python dev/serve.py            # serves site/ on http://localhost:8000

GitHub Pages answers range requests; Python's built-in http.server does
not, so the browser's gpxil reader needs this to work locally.
"""

from __future__ import annotations

import argparse
import os
import re
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

_RANGE = re.compile(r"bytes=(\d*)-(\d*)$")


class RangeHandler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map, ".js": "text/javascript", ".tif": "image/tiff", ".gz": "application/gzip"}

    def send_head(self):
        rng = self.headers.get("Range")
        path = self.translate_path(self.path)
        if not rng or os.path.isdir(path) or not os.path.isfile(path):
            return super().send_head()
        m = _RANGE.match(rng.strip())
        size = os.path.getsize(path)
        if not m:
            self.send_error(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            return None
        a, b = m.groups()
        if a == "":  # suffix range: the last b bytes
            start, end = max(0, size - int(b)), size - 1
        else:
            start, end = int(a), min(size - 1, int(b) if b else size - 1)
        if start >= size or start > end:
            self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return None
        fh = open(path, "rb")
        fh.seek(start)
        self.send_response(HTTPStatus.PARTIAL_CONTENT)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self._remaining = end - start + 1
        return fh

    def copyfile(self, source, outputfile):
        remaining = getattr(self, "_remaining", None)
        if remaining is None:
            return super().copyfile(source, outputfile)
        while remaining > 0:
            chunk = source.read(min(1 << 16, remaining))
            if not chunk:
                break
            outputfile.write(chunk)
            remaining -= len(chunk)
        self._remaining = None

    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, fmt, *args):  # quiet
        pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "site"))
    a = ap.parse_args()
    handler = partial(RangeHandler, directory=os.path.abspath(a.dir))
    with ThreadingHTTPServer(("127.0.0.1", a.port), handler) as httpd:
        print(f"serving {os.path.abspath(a.dir)} on http://localhost:{a.port}")
        httpd.serve_forever()


if __name__ == "__main__":
    main()
