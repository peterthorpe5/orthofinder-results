# Offline protein viewers

Runtime assets are already packaged; Node is not needed to use the app.
Maintainers can run `npm ci` then `npm run build` in this directory to rebuild
the pinned viewers. The build checks the pinned JBrowse checksum and records
asset checksums and third-party notices in the package data directory.

`npm test` validates the Alignment Viewer Sass export conversion. Generate
synthetic HTML with `python scripts/build_viewer_smoke.py --output /path/to/html`
from the repository root, then run `node frontend/browser-smoke.mjs /path/to/html`
with Playwright/Chromium installed. That check verifies real drawing and absence
of outbound HTTP requests. Auspice is a separately installed service.
