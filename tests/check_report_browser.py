#!/usr/bin/env python3
"""Exercise the generated standalone HTML in an offline Chromium context."""

from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def check(path):
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_context(offline=True).new_page()
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(path.as_uri())
        tp53 = page.locator("#summary tbody tr").nth(0)
        apc = page.locator("#summary tbody tr").nth(1)
        expect(tp53.locator("td").first).to_have_text("1 / 4 (25.0%)")
        expect(apc.locator("td").first).to_have_text("1 / 5 (20.0%)")

        cell = page.get_by_role("button", name="SYNTHETIC-02 · APC · GDC copy · Copy loss", exact=True)
        cell.focus()
        cell.press("Enter")
        expect(cell).to_have_attribute("aria-pressed", "true")
        expect(page.locator("#selection")).to_contain_text("SYNTHETIC-02 · APC")
        expect(page.locator("#details")).to_contain_text("mean -0.25")
        expect(page.locator("#details")).to_contain_text("Nonsense_Mutation")

        page.locator("#threshold").select_option("0.15")
        expect(tp53.locator("td").first).to_have_text("2 / 4 (50.0%)")
        expect(page.get_by_role("button", name="SYNTHETIC-01 · TP53 · GDC copy · Copy loss",
                                exact=True)).to_have_count(1)
        expect(page.get_by_role("button", name="SYNTHETIC-06 · TP53 · GDC copy · Copy gain",
                                exact=True)).to_have_count(1)

        page.locator("#threshold").select_option("0.25")
        expect(tp53.locator("td").first).to_have_text("0 / 4 (0.0%)")
        expect(apc.locator("td").first).to_have_text("1 / 5 (20.0%)")
        expect(page.get_by_role("button", name="SYNTHETIC-03 · TP53 · GDC copy · Copy loss",
                                exact=True)).to_have_count(1)
        expect(page.get_by_role("button", name="SYNTHETIC-04 · APC · GDC copy · Copy gain",
                                exact=True)).to_have_count(1)
        page.locator("#scope").select_option("matched")
        expect(apc.locator("td").first).to_have_text("1 / 4 (25.0%)")
        expect(tp53.locator("td").last).to_have_text("2 / 4 (50.0%)")
        assert not errors, errors


if __name__ == "__main__":
    with TemporaryDirectory() as directory:
        output = Path(directory) / "report.html"
        subprocess.run([sys.executable, "-B", str(ROOT / "examples/showcase/run_report.py"),
                        "--output", str(output)], check=True)
        check(output)
