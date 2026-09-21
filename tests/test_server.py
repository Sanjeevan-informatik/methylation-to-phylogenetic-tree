"""Local browser/API boundary and reproducible download integration."""

import base64
import io
import json
import threading
import urllib.error
import urllib.request
import zipfile
from http.server import ThreadingHTTPServer

import pytest

from methylphylo.server import Handler


@pytest.fixture
def service():
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}"
    server.shutdown()
    server.server_close()
    thread.join()


def test_local_ui_run_and_export(service):
    with urllib.request.urlopen(service) as response:
        assert b"Build a research run" in response.read()
    with urllib.request.urlopen(service + "/api/catalog") as response:
        assert len(json.load(response)["recipes"]) == 17
    text = "locus,A,B,C,D\nx,0.1,0.2,0.8,0.9\ny,0.2,0.3,0.7,0.8\n"
    payload = {"text": text, "data_type": "beta", "recipe": "beta-nj"}
    request = urllib.request.Request(
        service + "/api/infer",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "X-MethylPhylo-Request": "1"},
    )
    with urllib.request.urlopen(request) as response:
        result = json.load(response)
    assert result["report"]["data"]["sites"] == 2
    archive = zipfile.ZipFile(io.BytesIO(base64.b64decode(result["zip"])))
    assert archive.read("source-input.txt").decode() == text
    assert {"tree.nwk", "tree.svg", "report.json", "retained.tsv"} <= set(archive.namelist())


@pytest.mark.parametrize(
    "headers",
    [
        {},
        {"X-MethylPhylo-Request": "1", "Origin": "https://external.example"},
        {"X-MethylPhylo-Request": "1", "Host": "external.example"},
    ],
)
def test_cross_origin_or_missing_app_header_rejected(service, headers):
    request = urllib.request.Request(service + "/api/infer", data=b"{}", headers=headers)
    with pytest.raises(urllib.error.HTTPError) as error:
        urllib.request.urlopen(request)
    assert error.value.code == 403
