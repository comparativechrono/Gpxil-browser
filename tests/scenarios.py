"""More end-to-end scenarios with screenshots (deep links, zoom levels, wrap, search kinds, mobile).

    python tests/scenarios.py [--url http://localhost:8000/]
"""

import argparse
import os
import sys

from playwright.sync_api import sync_playwright

SCENARIOS = [
    # name, hash, viewport, actions
    ("10-deltaF508-letters", "#chr7:117559560-117559620", (1440, 900), []),
    ("11-alpha-satellite-171", "#chr1:123400001-123430000&wrap=171", (1440, 900), []),
    ("12-chr1-whole", "#chr1", (1440, 900), []),
    ("13-band-search", "#genome", (1440, 900), [("search", "17q21.31")]),
    ("14-variant-rs334", "#genome", (1440, 900), [("search", "rs334"), ("wait-hash", "#genome"), ("wait", 1500)]),
    ("15-alias-HER2", "#genome", (1440, 900), [("search", "HER2")]),
    ("16-mobile", "#chr17:7661779-7687538&gene=ENSG00000141510&mode=transcript", (390, 844), []),
    ("17-neon-cpg", "#chr17:7686000-7688200&pal=neon", (1440, 900), [("check", "cpg")]),
    ("18-tp53-expanded", "#chr17:7661779-7687538&gene=ENSG00000141510&mode=expanded", (1440, 900), [("wait", 4000)]),
    ("19-find-period", "#chr1:123400001-123430000", (1440, 900), [("click", "#find-period"), ("wait", 2500)]),
    ("20-live-variants", "#chr7:117559400-117559800", (1440, 900), [("check-live", "variation"), ("wait-live", "variation"), ("hover", (400, 300))]),
    ("21-live-regulatory", "#chr17:7660000-7700000", (1440, 900), [("check-live", "regulatory"), ("wait", 1500), ("shot", "loading"), ("wait-live", "regulatory")]),
    ("22-selection", "#chr7:117480000-117482000", (1440, 900), [("shiftdrag", ((200, 200), (600, 260))), ("wait", 800)]),
    ("23-karyo-hover", "#genome", (1440, 900), [("hover-karyo", "chr17"), ("wait", 300)]),
    ("24-karyo-chrM", "#genome", (1440, 900), [("wait", 1500), ("hover-karyo", "chrM"), ("wait", 300)]),
    ("25-regulatory-wide", "#chr17:7000000-8000000", (1440, 900), [("check-live", "regulatory"), ("wait-live", "regulatory")]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8000/")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "output"))
    ap.add_argument("--only")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for name, hash_, (w, h), actions in SCENARIOS:
            if a.only and a.only not in name:
                continue
            ctx = browser.new_context(viewport={"width": w, "height": h}, device_scale_factor=1, is_mobile=w < 500, has_touch=w < 500)
            ctx.add_init_script("localStorage.setItem('gpxil.skipIntro', 'true')")
            page = ctx.new_page()
            page.on("console", lambda m, n=name: errors.append(f"{n}: console.{m.type}: {m.text}") if m.type in ("error", "warning") else None)
            page.on("pageerror", lambda e, n=name: errors.append(f"{n}: pageerror: {e}"))
            page.goto(a.url + hash_)
            page.wait_for_selector("#app:not([hidden])", timeout=20000)
            page.wait_for_timeout(2500)
            for kind, arg in actions:
                if kind == "search":
                    page.fill("#search", arg)
                    page.wait_for_timeout(1200)
                    page.keyboard.press("Enter")
                    page.wait_for_timeout(3500)
                elif kind == "check":
                    page.evaluate(f"window.gpxil.view.setOptions({{{arg}: true}})")
                    page.wait_for_timeout(800)
                elif kind == "wait":
                    page.wait_for_timeout(arg)
                elif kind == "click":
                    page.evaluate("document.querySelector('.tabs button[data-tab=view]').click()")
                    page.click(arg)
                elif kind == "check-live":
                    page.evaluate("document.querySelector('.tabs button[data-tab=view]').click()")
                    page.check(f"[data-live={arg}]")
                elif kind == "hover":
                    box = page.query_selector("#view").bounding_box()
                    page.mouse.move(box["x"] + arg[0], box["y"] + arg[1])
                    page.wait_for_timeout(400)
                elif kind == "shiftdrag":
                    box = page.query_selector("#view").bounding_box()
                    (x0, y0), (x1, y1) = arg
                    page.keyboard.down("Shift")
                    page.mouse.move(box["x"] + x0, box["y"] + y0)
                    page.mouse.down()
                    page.mouse.move(box["x"] + x1, box["y"] + y1, steps=8)
                    page.mouse.up()
                    page.keyboard.up("Shift")
                elif kind == "wait-live":
                    try:
                        page.wait_for_function(f"window.gpxil.live.data['{arg}'] !== null", timeout=90000)
                    except Exception:  # Ensembl can be slow or down; the screenshot shows what the user would see
                        errors.append(f"{name}: no {arg} data from Ensembl within 90 s")
                    page.wait_for_timeout(600)
                elif kind == "wait-hash":
                    page.wait_for_function(f"location.hash !== '{arg}'", timeout=90000)
                elif kind == "shot":
                    page.screenshot(path=os.path.join(a.out, f"{name}-{arg}.png"))
                elif kind == "hover-karyo":
                    pt = page.evaluate(f"(() => {{ const b = window.gpxil.karyo.boxes.find(b => b.contig.name === '{arg}'); return [b.x + b.w/2, b.y + b.h*0.4]; }})()")
                    box = page.query_selector("#karyo").bounding_box()
                    page.mouse.move(box["x"] + pt[0], box["y"] + pt[1])
            page.screenshot(path=os.path.join(a.out, f"{name}.png"))
            print(name, page.evaluate("location.hash"))
            ctx.close()
        browser.close()
    for e in errors:
        print(e)
    return 1 if any("pageerror" in e or "console.error" in e for e in errors) else 0


if __name__ == "__main__":
    sys.exit(main())
