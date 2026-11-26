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
  await view("All results");

  // Single acquisition: Spectra. A scan redraws the spectrum, a match picks its scan
  await page.getByText("Single acquisition", { exact: true }).click(); await settle(3000);
  await view("Spectra");
  frames = insight();
  assert.equal(frames.length, 5, "Spectra shows scans, spectrum, chromatogram, peak map and matches");
  const SPECTRUM = "Spectrum m/z Intensity";
  const SCANS = "Scans Index";
  const MATCHES = "Spectrum matches";
  const spectrum = await text(await panel(SPECTRUM));
  await (await panel(SCANS)).locator(".tabulator-row").nth(2).click(); await settle(3000);
  const afterScan = await text(await panel(SPECTRUM));
  assert.notEqual(afterScan, spectrum, "spectrum follows the picked scan");
  await (await panel(MATCHES)).locator(".tabulator-row").nth(3).click(); await settle(3000);
  const afterMatch = await text(await panel(SPECTRUM));
  assert.notEqual(afterMatch, afterScan, "spectrum follows the picked match");
  const picked = await (await panel(SCANS)).locator(".tabulator-row.tabulator-selected").allInnerTexts();
  assert.equal(picked.length, 1, "the match's scan is picked in the scan table");
  await noException("Spectra");

  // The pick stays when the view changes
  await view("Identification quality");
  assert.equal(insight().length, 3, "Identification quality shows scores, errors and matches");
  await noException("Identification quality");
  await view("Spectra");
  assert.equal(await text(await panel(SPECTRUM)), afterMatch, "the picked spectrum stays across views");
  assert.deepEqual(await (await panel(SCANS)).locator(".tabulator-row.tabulator-selected").allInnerTexts(),
    picked, "the picked scan stays across views");

  assert.deepEqual(errors, [], "no script errors in the page");
  await browser.close();
  console.log("linked panels: ok");
})().catch((e) => { console.error(e); process.exit(1); });
