function toggleDetail(id) {
  const detail = document.getElementById(id);
  if (!detail) return;
  const isOpen = detail.style.display === "block";
  detail.style.display = isOpen ? "none" : "block";
  const card = detail.closest(".software-box");
  if (card) {
    card.classList.toggle("expanded", !isOpen);
  }
}

function toggleBox(id) {
  document.querySelectorAll(".abstract-box, .bibtex").forEach((box) => {
    if (box.id !== id) {
      box.classList.add("noshow");
    }
  });

  const box = document.getElementById(id);
  if (box) {
    box.classList.toggle("noshow");
  }
}

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll(".bibtex").forEach((bib) => {
    if (bib.querySelector(".copy-btn")) return;

    const btn = document.createElement("button");
    btn.className = "copy-btn";
    btn.textContent = "📋";
    btn.onclick = () => {
      const pre = bib.querySelector("pre");
      if (!pre) return;
      navigator.clipboard.writeText(pre.innerText).then(() => {
        btn.textContent = "✔️";
        setTimeout(() => {
          btn.textContent = "📋";
        }, 1500);
      });
    };

    bib.appendChild(btn);
  });
});

/* Keyword and co-author filters on the Publications page.
   Within a dropdown the choices are additive (any of them matches); across the
   two dropdowns they combine, so a publication must satisfy both. Year headings
   hide themselves when nothing is left underneath. */
document.addEventListener("DOMContentLoaded", () => {
  const filter = document.getElementById("pub-filter");
  if (!filter) return;

  const items = Array.from(document.querySelectorAll(".pub-item"));
  const status = document.getElementById("pub-filter-status");

  // Build the co-author options from the publications themselves, keeping only
  // names that appear on at least two.
  const MIN_JOINT = 2;
  const tally = new Map();
  items.forEach((item) => {
    (item.dataset.authors || "").split("|").filter(Boolean).forEach((name) => {
      tally.set(name, (tally.get(name) || 0) + 1);
    });
  });
  const authors = Array.from(tally.entries())
    .filter(([, n]) => n >= MIN_JOINT)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));

  const authorOptions = document.getElementById("pub-author-options");
  authors.forEach(([name, n]) => {
    const label = document.createElement("label");
    label.className = "pub-filter__option";
    const box = document.createElement("input");
    box.type = "checkbox";
    box.value = name;
    const nameEl = document.createElement("span");
    nameEl.className = "pub-filter__option-name";
    nameEl.textContent = name;
    const countEl = document.createElement("span");
    countEl.className = "pub-filter__count";
    countEl.textContent = n;
    label.append(box, nameEl, countEl);
    authorOptions.appendChild(label);
  });

  const groups = [
    {
      toggle: document.getElementById("pub-filter-toggle"),
      text: document.getElementById("pub-filter-toggle-text"),
      panel: document.getElementById("pub-filter-panel"),
      clear: document.getElementById("pub-filter-clear"),
      field: "keywords",
      empty: "all publications",
      unit: "keyword",
    },
    {
      toggle: document.getElementById("pub-author-toggle"),
      text: document.getElementById("pub-author-toggle-text"),
      panel: document.getElementById("pub-author-panel"),
      clear: document.getElementById("pub-author-clear"),
      field: "authors",
      empty: "anyone",
      unit: "co-author",
    },
  ];

  groups.forEach((g) => {
    g.boxes = Array.from(g.panel.querySelectorAll("input[type=checkbox]"));
  });

  const selectedOf = (g) => g.boxes.filter((b) => b.checked).map((b) => b.value);

  function setOpen(g, open) {
    g.panel.hidden = !open;
    g.toggle.setAttribute("aria-expanded", String(open));
    g.toggle.closest(".pub-filter__dropdown").classList.toggle("is-open", open);
    if (open) groups.filter((o) => o !== g).forEach((o) => setOpen(o, false));
  }

  function apply() {
    const picks = groups.map(selectedOf);
    let shown = 0;

    items.forEach((item) => {
      const ok = groups.every((g, i) => {
        if (picks[i].length === 0) return true;
        const hay = item.dataset[g.field] || "";
        return picks[i].some((v) => hay.includes("|" + v + "|"));
      });
      item.classList.toggle("is-hidden", !ok);
      if (ok) shown++;
    });

    document.querySelectorAll("[data-year-label]").forEach((label) => {
      let el = label.nextElementSibling;
      let keep = false;
      while (el && !el.hasAttribute("data-year-label")) {
        if (el.classList.contains("pub-item") && !el.classList.contains("is-hidden")) {
          keep = true;
          break;
        }
        el = el.nextElementSibling;
      }
      label.classList.toggle("is-hidden", !keep);
    });

    groups.forEach((g, i) => {
      const sel = picks[i];
      g.text.textContent =
        sel.length === 0 ? g.empty
          : sel.length === 1 ? sel[0]
          : sel.length + " " + g.unit + "s";
      g.toggle.closest(".pub-filter__dropdown").classList.toggle("has-selection", sel.length > 0);
    });

    const all = picks.flat();
    status.textContent = all.length === 0 ? ""
      : shown + (shown === 1 ? " publication" : " publications") + " matching " + all.join(", ");
  }

  groups.forEach((g) => {
    g.toggle.addEventListener("click", () => setOpen(g, g.panel.hidden));
    g.boxes.forEach((b) => b.addEventListener("change", apply));
    g.clear.addEventListener("click", () => {
      g.boxes.forEach((b) => (b.checked = false));
      apply();
    });
  });

  document.addEventListener("click", (e) => {
    if (!filter.contains(e.target)) groups.forEach((g) => setOpen(g, false));
  });

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      const open = groups.find((g) => !g.panel.hidden);
      if (open) {
        setOpen(open, false);
        open.toggle.focus();
      }
    }
  });

  apply();
});

/* Role filter on the Students page. Single-select: one level at a time, or
   "all". Year headings hide themselves when nothing is left underneath. */
document.addEventListener("DOMContentLoaded", () => {
  const filter = document.getElementById("stud-filter");
  if (!filter) return;

  const chips = Array.from(filter.querySelectorAll(".stud-filter__chip"));
  const items = Array.from(document.querySelectorAll(".stud-item"));

  function apply(role) {
    items.forEach((item) => {
      item.classList.toggle("is-hidden", role !== "" && item.dataset.role !== role);
    });

    document.querySelectorAll("[data-stud-year]").forEach((label) => {
      let el = label.nextElementSibling;
      let keep = false;
      while (el && !el.hasAttribute("data-stud-year")) {
        if (el.classList.contains("stud-item") && !el.classList.contains("is-hidden")) {
          keep = true;
          break;
        }
        el = el.nextElementSibling;
      }
      label.classList.toggle("is-hidden", !keep);
    });

    chips.forEach((c) => c.classList.toggle("is-active", c.dataset.role === role));
  }

  chips.forEach((chip) => {
    chip.addEventListener("click", () => apply(chip.dataset.role));
  });

  apply("");
});

/* Category filter on the Organized workshops page. Same shape as the Students
   role filter: one category at a time, or "all", with year headings hiding
   themselves when nothing is left underneath. */
document.addEventListener("DOMContentLoaded", () => {
  const filter = document.getElementById("event-filter");
  if (!filter) return;

  const chips = Array.from(filter.querySelectorAll(".event-filter__chip"));
  const cards = Array.from(document.querySelectorAll(".event-card[data-category]"));

  function apply(category) {
    cards.forEach((card) => {
      card.classList.toggle(
        "is-hidden",
        category !== "" && card.dataset.category !== category
      );
    });

    document.querySelectorAll("[data-event-year]").forEach((label) => {
      let el = label.nextElementSibling;
      let keep = false;
      while (el && !el.hasAttribute("data-event-year")) {
        if (el.classList.contains("event-card") && !el.classList.contains("is-hidden")) {
          keep = true;
          break;
        }
        el = el.nextElementSibling;
      }
      label.classList.toggle("is-hidden", !keep);
    });

    chips.forEach((c) => c.classList.toggle("is-active", c.dataset.category === category));
  }

  chips.forEach((chip) => {
    chip.addEventListener("click", () => apply(chip.dataset.category));
  });

  apply("");
});
