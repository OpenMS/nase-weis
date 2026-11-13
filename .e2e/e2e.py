"""Browser end-to-end test of NASEWEIS: every page, the example search and the decoy builder."""
import re
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8501"
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "e2e-output")
OUT.mkdir(parents=True, exist_ok=True)
FAIL = []
n = 0


def settle(page, timeout=180, quiet=2.0):
    start, idle = time.time(), None
    while time.time() - start < timeout:
        state = page.locator('[data-testid="stApp"]').first.get_attribute("data-test-script-state")
        if state == "notRunning":
            idle = idle or time.time()
            if time.time() - idle >= quiet:
                return
        else:
            idle = None
        page.wait_for_timeout(250)
    print(f"WARN: page still running after {timeout}s: {page.url}")


def shot(page, name):
    global n
    n += 1
    path = OUT / f"{n:02d}-{name}.png"
    page.screenshot(path=str(path), full_page=True)
    errs = page.locator('[data-testid="stException"]').all_inner_texts()
    for e in errs:
        FAIL.append(f"{name}: {e[:400]}")
    print(f"shot {path.name} exceptions={len(errs)}")


def goto(page, path, name):
    page.goto(f"{BASE}/{path}?workspace=5f0c8e2a-1d2b-4c3e-9f4a-7b6c5d4e3f21", wait_until="domcontentloaded")
    page.wait_for_selector('[data-testid="stApp"]', timeout=120000)
    settle(page)
    shot(page, name)


def click(page, text, exact=False):
    loc = page.get_by_role("button", name=re.compile(text) if not exact else text)
    if loc.count() == 0:
        FAIL.append(f"button not found: {text} on {page.url}")
        return False
    loc.first.click()
    settle(page)
    return True


def tab(page, text):
    loc = page.get_by_role("tab", name=text)
    if loc.count():
        loc.first.click()
        page.wait_for_timeout(1500)
        settle(page)
        return True
    FAIL.append(f"tab not found: {text} on {page.url}")
    return False


def wait_run(page, name, minutes=30):
    deadline = time.time() + minutes * 60
    k = 0
    while time.time() < deadline:
        body = page.locator("body").inner_text()
        if "completed successfully" in body:
            shot(page, f"{name}-finished")
            return True
        if re.search(r"Workflow (failed|cancelled)|ERROR", body) and "Stop Workflow" not in body:
            shot(page, f"{name}-failed")
            FAIL.append(f"{name}: run did not finish successfully")
            print(body[-3000:])
            return False
        k += 1
        if k % 12 == 1:
            shot(page, f"{name}-running")
        page.wait_for_timeout(5000)
    FAIL.append(f"{name}: run timed out")
    shot(page, f"{name}-timeout")
    return False


with sync_playwright() as p:
    for theme in ["light"]:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport={"width": 1440, "height": 900}, record_video_dir=str(OUT / "video"),
                                  record_video_size={"width": 1440, "height": 900}, color_scheme=theme)
        page = ctx.new_page()
        goto(page, "", "overview")
        goto(page, "topp_workflow_file_upload", "upload")
        click(page, "Load example dataset")
        shot(page, "upload-example-loaded")
        goto(page, "naseweis", "configure-identification")
        page.get_by_text("Advanced").first.click(); settle(page)
        shot(page, "configure-advanced")
        tab(page, "Quantification"); shot(page, "configure-quantification")
        goto(page, "topp_workflow_execution", "run")
        if click(page, "Start Workflow"):
            wait_run(page, "search")
        goto(page, "topp_workflow_results", "results-matches")
        # Experiment tabs, then one acquisition's own tabs under Sample details.
        for t in ["Sequence coverage", "Quantification", "Sample details", "Spectra", "Export"]:
            if tab(page, t):
                shot(page, "results-" + t.lower().replace(" ", "-"))
        goto(page, "naseweis_raw_viewer", "raw-viewer")
        goto(page, "decoy_database", "decoy-upload")
        inputs = page.locator('input[type="file"]')
        if inputs.count():
            inputs.first.set_input_files(sys.argv[3] if len(sys.argv) > 3 else "example-data/fasta/NME1_Target.fasta")
            page.wait_for_timeout(2000)
            click(page, "Add")
            shot(page, "decoy-uploaded")
        for step in ["Configure", "Run"]:
            page.locator('[data-testid="stRadio"]').get_by_text(step, exact=True).first.click(); settle(page)
            shot(page, f"decoy-{step.lower()}")
        if click(page, "Start Workflow"):
            wait_run(page, "decoy", minutes=10)
        page.locator('[data-testid="stRadio"]').get_by_text("Results", exact=True).first.click(); settle(page)
        shot(page, "decoy-results")
        ctx.close(); browser.close()

# dark theme pass over results
with sync_playwright() as p:
    browser = p.chromium.launch()
    ctx = browser.new_context(viewport={"width": 1440, "height": 900}, color_scheme="dark")
    page = ctx.new_page()
    goto(page, "topp_workflow_results", "dark-results")
    goto(page, "naseweis", "dark-configure")
    browser.close()

print("=" * 60)
if FAIL:
    print("E2E FAILURES:")
    for f in FAIL:
        print(" -", f)
else:
    print("E2E PASSED")
(OUT / "result.txt").write_text("\n".join(FAIL) or "PASSED")
sys.exit(1 if FAIL else 0)
