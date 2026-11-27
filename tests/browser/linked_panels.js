// Browser check of the linked panels of Explore results.
//
//   node tests/browser/linked_panels.js http://localhost:8501
//
// Needs Playwright with a Chromium; tests/test_naseweis_browser_links.py starts
// the app on the demo workspace and runs this script. Exits non-zero on the
// first panel that does not follow a pick.
const assert = require("node:assert/strict");
const { chromium } = require("playwright");

const url = process.argv[2] || "http://localhost:8501";

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1400, height: 1100 } });
  const errors = [];
  page.on("pageerror", (e) => errors.push(String(e)));
  const settle = async (ms = 1500) => {
    await page.waitForTimeout(ms);
    await page.waitForFunction(() => !document.querySelector('[data-testid="stStatusWidget"]'),
      null, { timeout: 90000 }).catch(() => {});
    await page.waitForTimeout(ms);
  };
  const insight = () => page.frames().filter((f) => /openms_insight/.test(f.url()));
  const viewBar = () => page.frames().find((f) => /view_bar/.test(f.url()));
  // Insight frames, found by what they draw (their order can change on a rerun)
  const panel = async (marker) => {
    for (const f of insight()) {
      const t = (await f.locator("body").innerText().catch(() => "")).replace(/\s+/g, " ");
      if (t.includes(marker)) return f;
    }
    throw new Error(`no panel shows "${marker}"`);
  };
  const text = async (frame) => (await frame.locator("body").innerText()).replace(/\s+/g, " ");
  const view = async (name) => { await viewBar().getByText(name, { exact: true }).click(); await settle(3000); };
  const noException = async (where) => {
    const n = await page.locator('[data-testid="stException"], [data-testid="stAlertContentError"]').count();
    assert.equal(n, 0, `${where}: the page shows an error`);
  };

  await page.goto(url); await settle();
  await page.getByRole("link", { name: /Explore results/ }).first().click(); await settle(2500);
  await noException("Explore results");

  // Experiment: Abundances. Picking an oligonucleotide redraws its profile
  await view("Abundances");
  let frames = insight();
  assert.equal(frames.length, 3, "Abundances shows table, profile and z-scores");
  const profile = await text(frames[1]);
  const abundanceRows = frames[0].locator(".tabulator-row");
  assert.ok(await abundanceRows.count() > 1, "abundance table has rows");
  await abundanceRows.nth(1).click(); await settle(3000);
  assert.notEqual(await text(insight()[1]), profile, "abundance profile follows the picked oligonucleotide");
  await noException("Abundances");
  await view("Identification");
  assert.ok(insight().length >= 1, "Identification shows the table of identified RNAs");
  await noException("Identification");
  await view("Summary");

  // Single acquisition: Spectra. A scan redraws the spectrum, a match picks its scan
  await page.getByText("Single acquisition", { exact: true }).click(); await settle(3000);
  await view("Spectra");
  frames = insight();
  assert.equal(frames.length, 5, "Spectra shows scans, spectrum, chromatogram, peak map and matches");
  const SPECTRUM = "Spectrum m/z Intensity";
  const SCANS = "Scans Index";
  const MATCHES = "Spectrum matches";
  // A spectrum is drawn on a canvas, and two scans can share their axis ticks,
  // so what the panel looks like (a screenshot) tells two spectra apart. Two
  // drawings of one spectrum can differ in a few anti-aliased pixels.
  const drawn = async () => {
    const f = await panel(SPECTRUM);
    await page.mouse.move(0, 0);  // no hover marks in the picture
    return (await f.locator("body").screenshot({ animations: "disabled" })).toString("base64");
  };
  const differing = (a, b) => page.evaluate(async ([a, b]) => {
    const pixels = async (png) => {
      const img = new Image(); img.src = "data:image/png;base64," + png; await img.decode();
      const c = document.createElement("canvas"); c.width = img.width; c.height = img.height;
      const g = c.getContext("2d"); g.drawImage(img, 0, 0);
      return g.getImageData(0, 0, c.width, c.height);
    };
    const [p, q] = [await pixels(a), await pixels(b)];
    if (p.width !== q.width || p.height !== q.height) return 1;
    let n = 0;
    for (let i = 0; i < p.data.length; i += 4) {
      const d = Math.abs(p.data[i] - q.data[i]) + Math.abs(p.data[i + 1] - q.data[i + 1])
        + Math.abs(p.data[i + 2] - q.data[i + 2]);
      if (d > 60) n++;
    }
    return n / (p.width * p.height);
  }, [a, b]);
  const same = async (a, b) => (await differing(a, b)) < 0.002;
  const pickUntilChanged = async (marker, from) => {
    for (let row = 1; row < 8; row++) {
      // A narrow table folds columns into the row; a click on its first cell picks it.
      await (await panel(marker)).locator(".tabulator-row").nth(row).locator(".tabulator-cell").first().click();
      // A slow machine redraws later: wait up to 15 s for the spectrum to change.
      for (let waited = 0; waited < 15000; waited += 500) {
        await settle(500);
        const now = await drawn();
        if (!(await same(now, from))) { await settle(1000); return await drawn(); }
      }
    }
    return from;
  };
  const spectrum = await drawn();
  const afterScan = await pickUntilChanged(SCANS, spectrum);
  assert.ok(!(await same(afterScan, spectrum)), "spectrum follows the picked scan");
  const afterMatch = await pickUntilChanged(MATCHES, afterScan);
  assert.ok(!(await same(afterMatch, afterScan)), "spectrum follows the picked match");
  const picked = await (await panel(SCANS)).locator(".tabulator-row.tabulator-selected").allInnerTexts();
  assert.equal(picked.length, 1, "the match's scan is picked in the scan table");
  await noException("Spectra");

  // The pick stays when the view changes
  await view("Search quality");
  assert.equal(insight().length, 3, "Search quality shows scores, errors and matches");
  await noException("Search quality");
  await view("Spectra");
  const back = await drawn();
  assert.ok(await same(back, afterMatch), `the picked spectrum stays across views (${await differing(back, afterMatch)})`);
  assert.deepEqual(await (await panel(SCANS)).locator(".tabulator-row.tabulator-selected").allInnerTexts(),
    picked, "the picked scan stays across views");

  assert.deepEqual(errors, [], "no script errors in the page");
  await browser.close();
  console.log("linked panels: ok");
})().catch((e) => { console.error(e); process.exit(1); });
