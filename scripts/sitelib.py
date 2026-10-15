"""Minimal static-site helpers for the project website.

No build dependency: every page is assembled from plain strings and written to
``docs/``, which GitHub Pages serves directly. Numbers are injected from
``results/tables`` at build time so the prose can never drift from the output.
"""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

AUTHOR = "Ashok Dwivedi"
AUTHOR_SITE = "https://ashokdwivedi.com"
REPO = "https://github.com/{owner}/{repo}"


@dataclass(frozen=True)
class Page:
    slug: str
    nav: str
    title: str
    subtitle: str


PAGES: tuple[Page, ...] = (
    Page("index.html", "Overview", "Mapping the tumour microenvironment in Xenium in situ data",
         "A reproducible, tested pipeline from raw 10x Xenium output to spatial "
         "statistics on the tumour\u2013immune interface."),
    Page("qc.html", "1 · Quality control", "Quality control: what to trust before anything else",
         "Per-cell, per-run and position-aware checks across two replicate sections, "
         "with explicit pass/fail criteria."),
    Page("celltypes.html", "2 · Cell typing", "Cell typing against a single-cell reference",
         "Harmony-integrated label transfer from a breast-cancer scRNA-seq atlas, "
         "cross-checked against canonical markers."),
    Page("spatial.html", "3 · Spatial structure", "Spatial neighbourhoods and the tumour\u2013immune interface",
         "Who touches whom, which recurrent niches the tissue is built from, and "
         "which populations are enriched at the tumour boundary."),
    Page("methods.html", "Methods", "Methods, reproducibility and limitations",
         "Every parameter, every threshold, and an honest account of what this "
         "analysis cannot tell you."),
    Page("about.html", "About", "About this project",
         "Why it was built, who built it, and how to get in touch."),
)


def _nav(current: str) -> str:
    items = []
    for p in PAGES:
        cls = ' class="current"' if p.slug == current else ""
        items.append(f'<a href="{p.slug}"{cls}>{p.nav}</a>')
    return "\n        ".join(items)


def render(page: Page, body: str, repo_url: str, *, description: str | None = None) -> str:
    """Wrap page ``body`` in the shared shell."""
    idx = [p.slug for p in PAGES].index(page.slug)
    prev_p = PAGES[idx - 1] if idx > 0 else None
    next_p = PAGES[idx + 1] if idx < len(PAGES) - 1 else None
    pager = ""
    if prev_p or next_p:
        left = (f'<a href="{prev_p.slug}"><span class="dir">Previous</span>'
                f"{prev_p.nav}</a>") if prev_p else "<span></span>"
        right = (f'<a href="{next_p.slug}" style="text-align:right"><span class="dir">Next</span>'
                 f"{next_p.nav}</a>") if next_p else "<span></span>"
        pager = f'<div class="pager">{left}{right}</div>'

    desc = escape(description or page.subtitle)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(page.title)} \u2014 Xenium breast TME</title>
<meta name="description" content="{desc}">
<meta name="author" content="{AUTHOR}">
<meta property="og:title" content="{escape(page.title)}">
<meta property="og:description" content="{desc}">
<meta property="og:type" content="article">
<link rel="stylesheet" href="assets/style.css">
</head>
<body>
<header class="masthead">
  <div class="masthead-inner">
    <a class="masthead-title" href="index.html">Xenium breast TME</a>
    <nav class="site-nav">
        {_nav(page.slug)}
        <a href="{repo_url}">GitHub \u2197</a>
    </nav>
  </div>
</header>
<main>
<p class="kicker">{escape(page.nav)}</p>
<h1>{escape(page.title)}</h1>
<p class="lede">{escape(page.subtitle)}</p>
{body}
{pager}
</main>
<footer class="site-footer">
  <div class="footer-inner">
    <p class="author-line">Built by <a href="{AUTHOR_SITE}">{AUTHOR}</a> \u2014
       <a href="{AUTHOR_SITE}">{AUTHOR_SITE.replace('https://', '')}</a></p>
    <p>Source, data-download script and unit tests:
       <a href="{repo_url}">{repo_url.replace('https://', '')}</a>. Released under the MIT licence.</p>
    <p>Data are public: 10x Genomics Xenium In Situ FFPE Human Breast Cancer (Janesick et al. 2023)
       and the breast-cancer scRNA-seq atlas of Wu et al. 2021 (GEO GSE176078). This site is an
       independent methods demonstration and is not affiliated with or endorsed by 10x Genomics.</p>
  </div>
</footer>
</body>
</html>
"""


def figure(src: str, caption: str, *, alt: str) -> str:
    return (f'<figure>\n<img src="{src}" alt="{escape(alt)}" loading="lazy">\n'
            f"<figcaption>{caption}</figcaption>\n</figure>")


def stats(pairs: list[tuple[str, str]]) -> str:
    cells = "".join(
        f'<div class="stat"><span class="v">{v}</span><span class="k">{escape(k)}</span></div>'
        for v, k in pairs
    )
    return f'<div class="stats">{cells}</div>'


def note(body: str, label: str = "Note", kind: str = "note") -> str:
    return f'<div class="{kind}"><span class="label">{escape(label)}</span>{body}</div>'


def caveat(body: str, label: str = "Limitation") -> str:
    return note(body, label, kind="caveat")


def table(headers: list[str], rows: list[list[str]], numeric_from: int = 1) -> str:
    head = "".join(
        f'<th class="num">{escape(h)}</th>' if i >= numeric_from else f"<th>{escape(h)}</th>"
        for i, h in enumerate(headers)
    )
    body = []
    for r in rows:
        tds = "".join(
            f'<td class="num">{c}</td>' if i >= numeric_from else f"<td>{c}</td>"
            for i, c in enumerate(r)
        )
        body.append(f"<tr>{tds}</tr>")
    return ('<div class="table-wrap"><table><thead><tr>' + head + "</tr></thead><tbody>"
            + "".join(body) + "</tbody></table></div>")


def cards(items: list[tuple[str, str, str, str]]) -> str:
    out = []
    for n, href, title, blurb in items:
        out.append(f'<a class="card" href="{href}"><span class="n">{escape(n)}</span>'
                   f"<h3>{escape(title)}</h3><p>{escape(blurb)}</p></a>")
    return f'<div class="cards">{"".join(out)}</div>'


def write(docs: Path, page: Page, body: str, repo_url: str) -> Path:
    out = docs / page.slug
    out.write_text(render(page, body, repo_url), encoding="utf-8")
    return out
