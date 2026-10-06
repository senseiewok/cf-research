"""Preview the CFTR model on this computer: serves web/ on 127.0.0.1 only, with a strict Content-Security-Policy.

Nothing else on the network can reach it. It sends the headers a strict host would (no inline script, no third-party request), so a mistake that such a policy
would block shows up here first. Every response
carries `Cache-Control: no-store`, so an edit appears on reload. Folders without an index are not listed, only files inside src/ are reachable,
and only GET and HEAD work. If the port is busy it tries the next few.

Usage: python tools/cftr-model/serve.py [--port 8080] [--no-browser]      (Ctrl+C or closing the window stops it)
Standard library only."""
import argparse
import functools
import http.server
import sys
import threading
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "web"
DEFAULT_PORT = 8080
CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; font-src 'self'; connect-src 'self'; "
       "base-uri 'none'; form-action 'none'; frame-ancestors 'none'; object-src 'none'")
# Python takes content types from the Windows registry, which can call .js "text/plain"; with nosniff the browser would then refuse the script.
TYPES = {".css": "text/css", ".js": "text/javascript", ".mjs": "text/javascript", ".json": "application/json", ".html": "text/html",
         ".webp": "image/webp", ".woff2": "font/woff2", ".txt": "text/plain", ".xml": "application/xml", ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon"}


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **TYPES}

    def end_headers(self):
        self.send_header("Content-Security-Policy", CSP)
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "strict-origin-when-cross-origin")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def list_directory(self, path):          # never list a folder; a missing index is a 404, as on the host
        self.send_error(404, "Not found")
        return None

    def send_error(self, code, message=None, explain=None):
        page = Path(self.directory) / "404.html"
        if code == 404 and page.is_file():
            body = page.read_bytes()
            self.send_response(404)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if self.command != "HEAD":
                self.wfile.write(body)
            return
        super().send_error(code, message, explain)

    def log_message(self, fmt, *args):
        sys.stderr.write("  %s\n" % (fmt % args))


def make_server(directory, port, tries=20):
    """A server for `directory` on 127.0.0.1, at `port` or the next free one (port 0 asks the system for any free port)."""
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"no such folder: {directory}")
    handler = functools.partial(Handler, directory=str(directory))
    last = None
    for p in ([0] if port == 0 else range(port, port + tries)):
        try:
            return http.server.ThreadingHTTPServer(("127.0.0.1", p), handler)
        except OSError as e:
            last = e
    raise OSError(f"no free port from {port} to {port + tries - 1}: {last}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Preview the site locally.")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args(argv)
    try:
        httpd = make_server(SRC, a.port)
    except (OSError, FileNotFoundError) as e:
        print(f"Could not start the preview: {e}")
        return 1
    url = f"http://127.0.0.1:{httpd.server_address[1]}/"
    print(f"CFTR model preview: {url}")
    print("It is served from web/ with a strict Content-Security-Policy, and only this computer can reach it.")
    print("Close this window (or press Ctrl+C) to stop the server.")
    sys.stdout.flush()
    if not a.no_browser:
        threading.Timer(0.6, webbrowser.open, args=(url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
