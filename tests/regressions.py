"""Regression checks in a real browser, for bugs found in review.

    python dev/serve.py &            # with the full data built
    python tests/regressions.py [--url http://localhost:8000/]

Each check prints PASS or FAIL; the exit code is the number of failures.
Ensembl answers are delayed or replaced where a check needs a race.
"""

from __future__ import annotations

import argparse
import json
import sys

from playwright.sync_api import sync_playwright

results: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, ok, detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'  (' + detail + ')' if detail else ''}", flush=True)


def new_page(browser, url, hash_="", *, skip_intro=True, viewport=(1440, 900), storage=None):
    ctx = browser.new_context(viewport={"width": viewport[0], "height": viewport[1]}, device_scale_factor=1)
    init = ""
    if skip_intro:
        init += "localStorage.setItem('gpxil.skipIntro', 'true');"
    for k, v in (storage or {}).items():
        init += f"localStorage.setItem({json.dumps(k)}, {json.dumps(json.dumps(v))});"
    if init:
        ctx.add_init_script(f"if (!sessionStorage.getItem('seeded')) {{ {init} sessionStorage.setItem('seeded', '1'); }}")
    page = ctx.new_page()
    errors: list[str] = []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.goto(url + hash_)
    page.wait_for_selector("#app:not([hidden])", timeout=30000)
    page.wait_for_timeout(1500)
    return ctx, page, errors


def status(page) -> str:
    return page.evaluate("document.querySelector('#st-zoom').textContent")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://localhost:8000/")
    a = ap.parse_args()
    url = a.url
    with sync_playwright() as p:
        browser = p.chromium.launch()

        # 1. a palette name every object has must not break the site, now or on the next visit
        ctx, page, errors = new_page(browser, url, "#genome&pal=toString")
        page.reload()
        page.wait_for_selector("#app:not([hidden])", timeout=30000)
        page.wait_for_timeout(800)
        check("pal=toString is ignored", not errors and page.evaluate("window.gpxil.palette") == "classic", "; ".join(errors))
        ctx.close()

        # 2 and 3. reloading a view keeps its zoom, in both wraps
        for start_hash in ["#chr7:117559400-117559800", "#chr17:7661779-7687538", "#chr1:123400001-123430000&wrap=171", "#chr1:1000000-1100000&wrap=8192", "#chr1"]:
            ctx, page, errors = new_page(browser, url, start_hash)
            first = (page.evaluate("location.hash"), status(page))
            seen = [first]
            for _ in range(4):
                page.reload()
                page.wait_for_selector("#app:not([hidden])", timeout=30000)
                page.wait_for_timeout(900)
                seen.append((page.evaluate("location.hash"), status(page)))
            zooms = {z for _, z in seen[1:]}
            whole = any(h.startswith("#chr1:1-") for h, _ in seen[1:]) and start_hash != "#chr1"
            check(f"reloads keep the zoom: {start_hash}", len(zooms) == 1 and not whole and not errors, f"{seen[0][1]} -> {seen[-1][1]}; {seen[-1][0]}")
            ctx.close()

        # 4. a slow Ensembl answer for one gene never lands on the next one
        ctx = browser.new_context(viewport={"width": 1440, "height": 900})
        ctx.add_init_script("localStorage.setItem('gpxil.skipIntro', 'true')")
        page = ctx.new_page()
        held = []

        def hold(route):
            if "ENSG00000141510" in route.request.url:  # TP53's transcripts: answer late
                held.append(route)
            else:
                route.continue_()

        page.route("https://rest.ensembl.org/lookup/id/**", hold)
        page.route("https://mygene.info/**", lambda r: r.fulfill(status=200, content_type="application/json", body=json.dumps({"summary": "SUMMARY FOR " + r.request.url.split("/")[-1].split("?")[0]})))
        page.goto(url + "#chr17:7661779-7687538&gene=ENSG00000141510&mode=transcript")
        page.wait_for_selector("#app:not([hidden])", timeout=30000)
        page.wait_for_timeout(1200)
        page.evaluate("window.gpxil.selectGene('ENSG00000012048', 'chr17', {mode: 'transcript'})")  # BRCA1
        page.wait_for_timeout(2500)
        for r in held:
            r.continue_()
        page.wait_for_timeout(4000)
        opts = page.evaluate("[...document.querySelectorAll('#gene-panel select[data-role=tx] option')].map(o => o.textContent)")
        title = page.evaluate("document.querySelector('#gene-panel h2').textContent")
        summary = page.evaluate("document.querySelector('#gene-panel .summary').textContent")
        check("late transcripts stay with their gene", "BRCA1" in title and opts and all("BRCA1" in o for o in opts) and "7157" not in summary, f"{title.strip()}: {opts[:2]}; {summary[:40]}")
        ctx.close()

        # 7. Enter goes where the box says, even before the suggestions catch up
        ctx, page, errors = new_page(browser, url, "#genome")
        page.fill("#search", "CFTR")
        page.wait_for_timeout(1500)
        page.fill("#search", "chr1:1,000,000-1,001,000")
        page.keyboard.press("Enter")
        page.wait_for_timeout(1500)
        check("Enter uses the text in the box", page.evaluate("location.hash").startswith("#chr1:"), page.evaluate("location.hash"))
        ctx.close()

        # 8. a gene highlighted on one chromosome is not written into another's link
        ctx, page, errors = new_page(browser, url, "#chr17:7661779-7687538&gene=ENSG00000141510")
        page.select_option("#contig-select", "chr1")
        page.wait_for_timeout(1200)
        check("no gene in another chromosome's link", "gene=" not in page.evaluate("location.hash"), page.evaluate("location.hash"))
        ctx.close()

        # 9. the palette survives a reload even with saved display options
        ctx, page, errors = new_page(browser, url, "#chr17:7686000-7688200")
        page.evaluate("document.querySelector('.tabs button[data-tab=view]').click()")
        page.check("[data-opt=cpg]")
        page.click("[data-pal=neon]")
        page.wait_for_timeout(500)
        page.goto(url + "#chr17:7686000-7688200")
        page.reload()
        page.wait_for_selector("#app:not([hidden])", timeout=30000)
        page.wait_for_timeout(1000)
        pal = page.evaluate("[window.gpxil.palette, window.gpxil.view.opts.palette, window.gpxil.view.opts.cpg]")
        check("palette and options both remembered", pal == ["neon", "neon", True], str(pal))
        ctx.close()

        # 10. a two-finger pinch is not a click
        ctx = browser.new_context(viewport={"width": 390, "height": 844}, is_mobile=True, has_touch=True)
        ctx.add_init_script("localStorage.setItem('gpxil.skipIntro', 'true')")
        page = ctx.new_page()
        page.goto(url + "#chr7:117480000-117670000")
        page.wait_for_selector("#app:not([hidden])", timeout=30000)
        page.wait_for_timeout(2000)
        box = page.query_selector("#view").bounding_box()
        cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
        cdp = ctx.new_cdp_session(page)

        def touch(kind, pts):
            cdp.send("Input.dispatchTouchEvent", {"type": kind, "touchPoints": [{"x": x, "y": y, "id": i} for i, (x, y) in enumerate(pts)]})

        touch("touchStart", [(cx, cy)])
        touch("touchStart", [(cx, cy), (cx + 10, cy + 10)])
        for k in range(1, 8):
            touch("touchMove", [(cx - 6 * k, cy - 6 * k), (cx + 10 + 6 * k, cy + 10 + 6 * k)])
        touch("touchEnd", [(cx - 42, cy - 42)])
        touch("touchEnd", [])
        page.wait_for_timeout(1200)
        check("pinch does not highlight a gene", page.evaluate("window.gpxil.highlight") is None)
        ctx.close()

        # 11. one-row rasters in the gene panel are visible, and hover maps to the right cell
        ctx, page, errors = new_page(browser, url, "#genome")
        gid = page.evaluate("(async () => { const r = await window.gpxil.genes.resolve('MIR6776'); return r && [r.row[1], r.row[2]]; })()")
        if gid:
            page.evaluate(f"window.gpxil.selectGene('{gid[0]}', '{gid[1]}', {{mode: 'transcript'}})")
            page.wait_for_timeout(3000)
            h = page.evaluate("(() => { const c = document.querySelector('#gene-panel canvas.mrna'); return c && [c.clientHeight, c.height]; })()")
            check("one-row mRNA raster is visible", bool(h) and h[0] >= h[1] and h[0] > 0, str(h))
        else:
            check("one-row mRNA raster is visible", False, "MIR6776 not found")
        ctx.close()

        # 12. scrolling a zoomed-out view (bases averaged per pixel) is quick after the first frame
        ctx, page, errors = new_page(browser, url, "#chr7:110000000-122000000")
        page.wait_for_timeout(3000)
        t = page.evaluate(
            """(async () => {
              const v = window.gpxil.view;
              const times = [];
              for (let i = 0; i < 8; i++) {
                v.scrollRows(5);
                const t0 = performance.now();
                v.buildFrame(v.layout());
                times.push(performance.now() - t0);
              }
              return [v.layout().bpp, Math.max(...times.slice(1))];
            })()"""
        )
        check("zoomed-out scroll frames are fast", t[0] > 1 and t[1] < 60, f"{t[0]} bases per pixel, slowest frame {t[1]:.0f} ms")
        ctx.close()

        # 13. fixed wrap widths are capped
        ctx, page, errors = new_page(browser, url, "#chr1:1000000-1100000&wrap=100000")
        check("fixed wrap is capped", page.evaluate("window.gpxil.view.state.width") <= 16384, str(page.evaluate("window.gpxil.view.state.width")))
        ctx.close()

        # 15. Space presses a focused button; Escape clears a linked point
        ctx, page, errors = new_page(browser, url, "#chr7:117559590")
        before = (status(page), page.evaluate("window.gpxil.view.state.scrollRow"))
        page.focus("#zoom-out")
        page.keyboard.press("Space")
        page.wait_for_timeout(600)
        after = (status(page), page.evaluate("window.gpxil.view.state.scrollRow"))
        check("Space presses the focused button", before[0] != after[0], f"{before[0]} -> {after[0]}")
        page.focus("#view")
        page.keyboard.press("Escape")
        page.wait_for_timeout(300)
        check("Escape clears a linked point", page.evaluate("window.gpxil.view.selection") is None)
        ctx.close()

        # 16. one bookmark click is one navigation
        ctx, page, errors = new_page(browser, url, "#chr17:7661779-7687538", storage={"gpxil.bookmarks": [{"label": "x", "hash": "#chr7:117559400-117559800", "when": "2026-10-27"}]})
        page.evaluate("window.gpxil.renderBookmarks(); window.__routes = 0; const r = window.gpxil.route.bind(window.gpxil); window.gpxil.route = (...a) => { window.__routes++; return r(...a); }; 0")  # (a trailing function would be called)
        page.evaluate("document.querySelector('.tabs button[data-tab=saved]').click()")
        page.click("#bookmarks a[data-bm='0']")
        page.wait_for_timeout(1200)
        check("bookmark click routes once", page.evaluate("window.__routes") == 1, f"{page.evaluate('window.__routes')} routes")
        ctx.close()

        # 17. "Zoom to gene" from the karyotype opens the gene's chromosome
        ctx, page, errors = new_page(browser, url, "#chr17:7661779-7687538&gene=ENSG00000141510")
        page.keyboard.press("k")
        page.wait_for_timeout(800)
        page.evaluate("document.querySelector('.tabs button[data-tab=gene]').click()")
        page.click("#gene-panel button[data-act=zoom]")
        page.wait_for_timeout(1200)
        st = page.evaluate("[window.gpxil.screen, document.querySelector('#contig-select').value]")
        check("Zoom to gene leaves the karyotype", st == ["contig", "chr17"], str(st))
        ctx.close()

        # Expanded shows every transcript as lanes, from the site's own data; clicking a lane opens that transcript
        ctx, page, errors = new_page(browser, url, "#chr17:7661779-7687538&gene=ENSG00000141510&mode=expanded")
        page.wait_for_timeout(1500)
        n = page.evaluate("window.gpxil.lanes.active && window.gpxil.lanes.txs.length")
        box = page.query_selector("#lanes").bounding_box()
        y = page.evaluate("window.gpxil.lanes.laneTop(3) + window.gpxil.lanes.laneH() / 2")
        picked = page.evaluate("window.gpxil.lanes.txs[2].id")
        page.mouse.click(box["x"] + box["width"] * 0.6, box["y"] + y)
        page.wait_for_timeout(1000)
        st = page.evaluate("[window.gpxil.panel.mode, window.gpxil.panel.tx.id, window.gpxil.lanes.active, window.gpxil.view.paused]")
        check("Expanded shows every transcript; a click opens one", bool(n) and n >= 30 and st == ["transcript", picked, False, False] and not errors, f"{n} lanes; {st}")
        ctx.close()

        # selenocysteine: GPX1's Sec49 is U, not a stop
        ctx, page, errors = new_page(browser, url, "#genome")
        r = page.evaluate("(async () => { const r = await window.gpxil.genes.resolve('GPX1'); return [r.row[1], r.row[2]]; })()")
        page.evaluate(f"window.gpxil.selectGene('{r[0]}', '{r[1]}', {{mode: 'transcript'}})")
        page.wait_for_timeout(3500)
        prot = page.evaluate("window.gpxil.panel.protein || ''")
        check("GPX1 reads selenocysteine", "U" in prot and "*" not in prot[:-1], f"U at {prot.find('U') + 1}, length {len(prot)}")
        ctx.close()

        browser.close()
    failed = sum(1 for _, ok, _ in results if not ok)
    print(f"{len(results) - failed} passed, {failed} failed")
    return failed


if __name__ == "__main__":
    sys.exit(main())
