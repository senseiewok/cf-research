"""The preview/test server must survive a burst of connections: a deep listen queue and daemon threads (the default backlog of 5 refuses connections on Windows).

Run from this folder:   python -m unittest test_serve -v"""
import socket
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import serve  # noqa: E402


class ServerBacklogTests(unittest.TestCase):
    def test_class_has_a_deep_listen_queue_and_daemon_threads(self):
        self.assertGreaterEqual(serve.Server.request_queue_size, 64)
        self.assertTrue(serve.Server.daemon_threads)

    def test_make_server_returns_that_class_and_listens_on_loopback_only(self):
        with tempfile.TemporaryDirectory() as d:
            httpd = serve.make_server(d, 0)
            try:
                self.assertIsInstance(httpd, serve.Server)
                self.assertEqual(httpd.server_address[0], "127.0.0.1")
                self.assertGreaterEqual(httpd.request_queue_size, 64)
            finally:
                httpd.server_close()


if __name__ == "__main__":
    unittest.main()
