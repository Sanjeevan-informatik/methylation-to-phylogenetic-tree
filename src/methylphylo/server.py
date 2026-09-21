"""Loopback-only UI for local research runs. Not a public hosting backend."""

from __future__ import annotations

import base64
import io
import json
import tempfile
import threading
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .datasets import load_dataset, prepare
from .inference import catalog, iqtree_binary
from .research import run_analysis

WEB = Path(__file__).parent / "web"
MAX_BODY = 12 * 1024 * 1024
RUN_LOCK = threading.Lock()


class Handler(BaseHTTPRequestHandler):
    def send(self, status, body, content_type="application/json"):
        body = body.encode() if isinstance(body, str) else body
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")
        self.send_header(
            "Content-Security-Policy",
            "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' blob:; object-src 'none'; frame-ancestors 'none'",
        )
        self.end_headers()
        self.wfile.write(body)

    def local_request(self):
        port = self.server.server_port
        hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}
        origin = self.headers.get("Origin")
        return self.headers.get("Host") in hosts and (
            origin is None or origin in {f"http://{host}" for host in hosts}
        )

    def do_GET(self):
        if not self.local_request():
            self.send(403, '{"error":"Local requests only"}')
            return
        if self.path == "/api/catalog":
            try:
                engine = iqtree_binary()
            except ValueError:
                engine = None
            self.send(200, json.dumps({"recipes": catalog(), "iqtree_available": bool(engine)}))
            return
        files = {
            "/": ("index.html", "text/html; charset=utf-8"),
            "/app.js": ("app.js", "text/javascript; charset=utf-8"),
            "/style.css": ("style.css", "text/css; charset=utf-8"),
        }
        if self.path not in files:
            self.send(404, '{"error":"Not found"}')
            return
        name, mime = files[self.path]
        self.send(200, (WEB / name).read_bytes(), mime)

    def do_POST(self):
        if not self.local_request() or self.headers.get("X-MethylPhylo-Request") != "1":
            self.send(403, '{"error":"Local application requests only"}')
            return
        if self.path != "/api/infer":
            self.send(404, '{"error":"Not found"}')
            return
        if not RUN_LOCK.acquire(blocking=False):
            self.send(409, '{"error":"An analysis is already running"}')
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_BODY:
                raise ValueError("Upload must be smaller than 12 MB; use the CLI for larger data")
            payload = json.loads(self.rfile.read(length))
            if not isinstance(payload, dict) or not isinstance(payload.get("text"), str):
                raise TypeError("Expected a text input and analysis settings")
            replicates = int(payload.get("bootstrap", 0))
            if not 0 <= replicates <= 25:
                raise ValueError(
                    "Browser bootstrap limit is 25; use the CLI for research-scale runs"
                )
            with tempfile.TemporaryDirectory(prefix="methylphylo-") as directory:
                root = Path(directory)
                path = root / "upload.txt"
                path.write_text(payload["text"])
                data = load_dataset(
                    path,
                    payload["data_type"],
                    genotype_policy=payload.get("genotype_policy", "homozygous"),
                )
                if len(data.names) > 150 or len(data.loci) > 100000:
                    raise ValueError("Browser limit: 150 samples and 100,000 sites; use the CLI")
                data = prepare(
                    data,
                    missing=payload.get("missing", "complete"),
                    sites=payload.get("sites", "all"),
                )
                output = root / "results"
                report = run_analysis(
                    data,
                    payload["recipe"],
                    output,
                    seed=int(payload.get("seed", 1)),
                    timeout=120,
                    replicates=replicates,
                )
                (output / "source-input.txt").write_text(payload["text"])
                bundle = io.BytesIO()
                with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
                    for item in output.rglob("*"):
                        if item.is_file() and item.suffix in {
                            ".txt",
                            ".json",
                            ".nwk",
                            ".svg",
                            ".tsv",
                            ".fasta",
                            ".log",
                            ".iqtree",
                        }:
                            archive.write(item, str(item.relative_to(output)))
                self.send(
                    200,
                    json.dumps(
                        {
                            "report": report,
                            "svg": (output / "tree.svg").read_text(),
                            "zip": base64.b64encode(bundle.getvalue()).decode(),
                        }
                    ),
                )
        except (ValueError, OSError, KeyError, TypeError) as exc:
            self.send(400, json.dumps({"error": str(exc)}))
        finally:
            RUN_LOCK.release()


def serve(port=8765):
    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print(f"MethylPhylo research workbench: http://127.0.0.1:{server.server_port}", flush=True)
    print(
        "Files are processed on this machine. Download results before closing the page.", flush=True
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
