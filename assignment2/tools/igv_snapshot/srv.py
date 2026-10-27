import http.server, os, re, sys
class H(http.server.SimpleHTTPRequestHandler):
    def send_head(self):
        path = self.translate_path(self.path)
        rng = self.headers.get("Range")
        if not rng or not os.path.isfile(path):
            return super().send_head()
        size = os.path.getsize(path)
        a, b = re.match(r"bytes=(\d*)-(\d*)", rng).groups()
        a = int(a or 0); b = min(int(b) if b else size - 1, size - 1)
        f = open(path, "rb"); f.seek(a)
        self.send_response(206)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Range", f"bytes {a}-{b}/{size}")
        self.send_header("Content-Length", str(b - a + 1))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self._left = b - a + 1
        return f
    def copyfile(self, src, dst):
        if hasattr(self, "_left"):
            dst.write(src.read(self._left))
        else:
            super().copyfile(src, dst)
    def log_message(self, *a): pass
http.server.ThreadingHTTPServer(("127.0.0.1", 8765), lambda *a, **k: H(*a, directory=sys.argv[1], **k)).serve_forever()
