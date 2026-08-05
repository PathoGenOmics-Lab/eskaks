// MathJax 3 config for pymdownx.arithmatex (generic mode) + Material instant nav.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

// Re-typeset after Material's instant navigation. Guarded so that if this fires
// before the MathJax library has loaded, it no-ops instead of throwing — an
// uncaught throw here would break Material's shared document$ chain (which also
// drives mermaid rendering).
document$.subscribe(() => {
  if (typeof MathJax !== "undefined" && MathJax.typesetPromise) {
    MathJax.startup?.document?.clear?.();
    MathJax.texReset?.();
    MathJax.typesetPromise?.();
  }
});
