"""End-to-end smoke test: drive the browser in headless Chromium and take screenshots.

    python dev/serve.py &        # serve site/ (with site/data built)
    python tests/smoke.py [--url http://localhost:8000/] [--out tests/output]
"""

import argparse
import os
import sys

from playwright.sync_api import sync_playwright


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8000/")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "output"))
    ap.add_argument("--width", type=int, default=1440)
    ap.add_argument("--height", type=int, default=900)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": a.width, "height": a.height})
        page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}") if m.type in ("error", "warning") else None)
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))

        def shot(name):
            page.screenshot(path=os.path.join(a.out, f"{name}.png"))
            print("screenshot", name)

        page.goto(a.url)
        page.wait_for_timeout(2500)
        shot("01-landing")
        page.click("#enter")
        page.wait_for_selector("#app:not([hidden])", timeout=20000)
        page.wait_for_timeout(1800)
        shot("02-karyotype")
        # open chromosome 7 by clicking near its middle
        box = page.evaluate("(() => { const b = window.gpxil.karyo.boxes.find(b => b.contig.name === 'chr7'); return [b.x + b.w/2, b.y + b.h*0.75]; })()")
        canvas = page.query_selector("#karyo").bounding_box()
        page.mouse.click(canvas["x"] + box[0], canvas["y"] + box[1])
        page.wait_for_timeout(3000)
        shot("03-chr7-region")
        # search a gene
        page.fill("#search", "CFTR")
        page.wait_for_timeout(1500)
        shot("04-suggest")
        page.keyboard.press("Enter")
        page.wait_for_timeout(4000)
        shot("05-cftr-gene")
        page.click(".modes button[data-mode=transcript]")
        page.wait_for_timeout(2500)
        shot("06-cftr-transcript")
        page.click(".modes button[data-mode=expanded]")
        page.wait_for_timeout(4000)
        shot("07-cftr-expanded")
        print("hash:", page.evaluate("location.hash"))
        browser.close()
    for e in errors:
        print(e)
    return 1 if any(e.startswith("pageerror") or "console.error" in e for e in errors) else 0


if __name__ == "__main__":
    sys.exit(main())
