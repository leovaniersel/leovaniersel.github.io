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

/* Keyword filter on the Publications page.
   Chips are additive: selecting several shows publications matching any of
   them. Year headings hide themselves when nothing is left underneath. */
document.addEventListener("DOMContentLoaded", () => {
  const filter = document.getElementById("pub-filter");
  if (!filter) return;

  const chips = Array.from(filter.querySelectorAll(".pub-filter__chip"));
  const allChip = chips.find((c) => c.dataset.keyword === "");
  const status = document.getElementById("pub-filter-status");
  const items = Array.from(document.querySelectorAll(".pub-item"));
  const selected = new Set();

  function apply() {
    let shown = 0;
    items.forEach((item) => {
      const kws = item.dataset.keywords || "";
      const match =
        selected.size === 0 ||
        Array.from(selected).some((kw) => kws.includes("|" + kw + "|"));
      item.classList.toggle("is-hidden", !match);
      if (match) shown++;
    });

    // A heading stays only if some visible publication follows it before the next heading.
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

    allChip.classList.toggle("is-active", selected.size === 0);
    if (selected.size === 0) {
      status.textContent = "";
    } else {
      status.textContent =
        shown + (shown === 1 ? " publication" : " publications") +
        " matching " + Array.from(selected).join(", ");
    }
  }

  chips.forEach((chip) => {
    chip.addEventListener("click", () => {
      const kw = chip.dataset.keyword;
      if (kw === "") {
        selected.clear();
        chips.forEach((c) => c.classList.remove("is-active"));
      } else if (selected.has(kw)) {
        selected.delete(kw);
        chip.classList.remove("is-active");
      } else {
        selected.add(kw);
        chip.classList.add("is-active");
      }
      apply();
    });
  });

  apply();
});
