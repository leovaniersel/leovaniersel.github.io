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
   A dropdown of checkboxes rather than a row of chips, because the keyword
   list outgrew the width of the page. Selections are additive: several
   keywords show publications matching any of them. Year headings hide
   themselves when nothing is left underneath. */
document.addEventListener("DOMContentLoaded", () => {
  const filter = document.getElementById("pub-filter");
  if (!filter) return;

  const toggle = document.getElementById("pub-filter-toggle");
  const toggleText = document.getElementById("pub-filter-toggle-text");
  const panel = document.getElementById("pub-filter-panel");
  const clearBtn = document.getElementById("pub-filter-clear");
  const status = document.getElementById("pub-filter-status");
  const boxes = Array.from(panel.querySelectorAll("input[type=checkbox]"));
  const items = Array.from(document.querySelectorAll(".pub-item"));

  function selectedKeywords() {
    return boxes.filter((b) => b.checked).map((b) => b.value);
  }

  function setOpen(open) {
    panel.hidden = !open;
    toggle.setAttribute("aria-expanded", String(open));
    filter.classList.toggle("is-open", open);
  }

  function apply() {
    const selected = selectedKeywords();
    let shown = 0;

    items.forEach((item) => {
      const kws = item.dataset.keywords || "";
      const match =
        selected.length === 0 ||
        selected.some((kw) => kws.includes("|" + kw + "|"));
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

    if (selected.length === 0) {
      toggleText.textContent = "all publications";
      status.textContent = "";
    } else {
      toggleText.textContent =
        selected.length === 1 ? selected[0] : selected.length + " keywords";
      status.textContent =
        shown + (shown === 1 ? " publication" : " publications") +
        " matching " + selected.join(", ");
    }
    filter.classList.toggle("has-selection", selected.length > 0);
  }

  toggle.addEventListener("click", () => setOpen(panel.hidden));
  boxes.forEach((b) => b.addEventListener("change", apply));

  clearBtn.addEventListener("click", () => {
    boxes.forEach((b) => (b.checked = false));
    apply();
  });

  document.addEventListener("click", (e) => {
    if (!filter.contains(e.target)) setOpen(false);
  });

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && !panel.hidden) {
      setOpen(false);
      toggle.focus();
    }
  });

  apply();
});
