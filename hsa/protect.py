"""Password-protect the static dashboard: the whole page is AES-256-GCM encrypted with a key derived from the
password (PBKDF2-SHA256), and a small login page decrypts it in the browser with WebCrypto.

The password is read from .dashboard_password (gitignored) and is never written into the repo.
"""
import base64
import json
import os

from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

from .config import ROOT

PASSWORD_FILE = ROOT / ".dashboard_password"
ITERATIONS = 600_000


def encrypt(html: str, password: str) -> dict:
    salt, iv = os.urandom(16), os.urandom(12)
    key = PBKDF2HMAC(algorithm=hashes.SHA256(), length=32, salt=salt, iterations=ITERATIONS).derive(password.encode())
    ct = AESGCM(key).encrypt(iv, html.encode(), None)        # ciphertext || 16-byte tag, as WebCrypto expects
    b64 = lambda b: base64.b64encode(b).decode()
    return {"salt": b64(salt), "iv": b64(iv), "ct": b64(ct), "iter": ITERATIONS}


LOGIN = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<link rel="icon" href="__FAVICON__">
<title>Human Spatial Atlas</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Schibsted+Grotesk:wght@500;700&family=Source+Sans+3:wght@400;600&display=swap">
<style>
*,*::before,*::after{box-sizing:border-box}
:root{--bg:#fbfbfc;--panel:#fff;--fg:#16181d;--muted:#5b6170;--line:#e4e6ea;--accent:#3454a3;--eosin:#c24d7a;--bad:#b2412f;color-scheme:light;
  --r-sm:6px;--r:10px;--r-lg:14px;
  --f-display:"Schibsted Grotesk","Helvetica Neue",Arial,sans-serif;--f-body:"Source Sans 3","Segoe UI",Roboto,sans-serif}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0f1115;--panel:#181b21;--fg:#e6e8ec;--muted:#9aa1ad;--line:#2a2e37;--accent:#8fa8dd;--eosin:#ec8fb2;--bad:#ef8a77;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#0f1115;--panel:#181b21;--fg:#e6e8ec;--muted:#9aa1ad;--line:#2a2e37;--accent:#8fa8dd;--eosin:#ec8fb2;--bad:#ef8a77;color-scheme:dark}
html,body{height:100%}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 var(--f-body);display:grid;place-items:center;padding-inline:16px}
form{width:100%;max-width:360px;background:var(--panel);border:1px solid var(--line);border-radius:var(--r-lg);padding:28px 24px;display:grid;gap:14px;justify-items:center;text-align:center}
.logo{width:72px;height:72px}
h1{font:700 1.4rem/1.15 var(--f-display);margin:0}
h1 span{color:var(--accent)}
p{margin:0;color:var(--muted);font-size:.9rem}
label{justify-self:stretch;text-align:left;font-size:.82rem;color:var(--muted);display:grid;gap:6px}
input{font:inherit;padding:9px 11px;border-radius:var(--r);border:1px solid var(--line);background:var(--bg);color:var(--fg);width:100%}
button{justify-self:stretch;font:600 .95rem var(--f-body);padding:10px;border-radius:var(--r);border:0;background:var(--accent);color:var(--panel);cursor:pointer}
input:focus-visible,button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.err{color:var(--bad);font-size:.85rem;min-height:1.2em}
.build{color:var(--muted);font-size:.72rem;margin-top:1.1rem;font-variant-numeric:tabular-nums;letter-spacing:.01em}
</style></head>
<body>
<form id="f" autocomplete="off">
  __LOGO__
  <h1>Human Spatial Atlas <span>HSA</span></h1>
  <p>This dashboard is private. Enter the password to view it.</p>
  <label for="pw">Password<input id="pw" type="password" required autofocus></label>
  <button type="submit" id="go">View dashboard</button>
  <div class="err" id="err" role="alert"></div>
  <div class="build" title="If this is older than you expect, your browser is serving a cached copy -- hard-reload">__BUILD__</div>
</form>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
(function () {
  const P = JSON.parse(document.getElementById("payload").textContent);
  const b = (s) => Uint8Array.from(atob(s), (c) => c.charCodeAt(0));
  async function open(pw) {
    const base = await crypto.subtle.importKey("raw", new TextEncoder().encode(pw), "PBKDF2", false, ["deriveKey"]);
    const key = await crypto.subtle.deriveKey({ name: "PBKDF2", salt: b(P.salt), iterations: P.iter, hash: "SHA-256" },
      base, { name: "AES-GCM", length: 256 }, false, ["decrypt"]);
    const plain = await crypto.subtle.decrypt({ name: "AES-GCM", iv: b(P.iv) }, key, b(P.ct));
    return new TextDecoder().decode(plain);
  }
  function show(html) { document.open(); document.write(html); document.close(); }
  let saved = null;
  try { saved = sessionStorage.getItem("hsa-pw"); } catch (e) {}
  if (saved) open(saved).then(show).catch(() => { try { sessionStorage.removeItem("hsa-pw"); } catch (e) {} });
  document.getElementById("f").addEventListener("submit", async (e) => {
    e.preventDefault();
    const pw = document.getElementById("pw").value, err = document.getElementById("err"), go = document.getElementById("go");
    go.disabled = true; go.textContent = "Unlocking"; err.textContent = "";
    try {
      const html = await open(pw);
      try { sessionStorage.setItem("hsa-pw", pw); } catch (e2) {}
      show(html);
    } catch (e2) {
      err.textContent = "That password is not correct. Try again.";
      go.disabled = false; go.textContent = "View dashboard";
    }
  });
})();
</script>
</body></html>
"""


def _build_stamp() -> str:
    """Date and commit of the build, shown on the login page.

    The dashboard is one encrypted blob, so a stale browser cache is invisible:
    the page looks identical and a plaintext search for new markup finds nothing
    because the content is ciphertext. This puts the build where it can be read
    before logging in.
    """
    import datetime as _dt
    import subprocess as _sp
    when = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        sha = _sp.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                      capture_output=True, text=True, timeout=5)
        rev = sha.stdout.strip() if sha.returncode == 0 else ""
    except Exception:                                           # noqa: BLE001
        rev = ""
    return f"build {when}" + (f" \u00b7 {rev}" if rev else "")


def protect(full_html: str, favicon: str, logo_svg: str) -> str | None:
    if not PASSWORD_FILE.exists():
        return None
    pw = PASSWORD_FILE.read_text().strip()
    payload = json.dumps(encrypt(full_html, pw))
    return (LOGIN.replace("__PAYLOAD__", payload).replace("__FAVICON__", favicon)
            .replace("__LOGO__", logo_svg).replace("__BUILD__", _build_stamp()))
