// Documenter uses RequireJS, whose global define() causes Mermaid's UMD build
// to register as an AMD module instead of setting window.mermaid — so
// auto-init never fires. Injecting a type="module" script bypasses RequireJS.
// We also call mermaid.run() explicitly rather than relying on startOnLoad,
// because DOMContentLoaded has already fired by the time a dynamically-injected
// module script executes.
(function () {
  var s = document.createElement('script');
  s.type = 'module';
  s.textContent = [
    'import mermaid from "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs";',
    'mermaid.initialize({ startOnLoad: false, theme: "default" });',
    'mermaid.run({ querySelector: ".mermaid" });',
  ].join('\n');
  document.head.appendChild(s);
})();
