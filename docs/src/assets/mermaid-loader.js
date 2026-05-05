// Documenter uses RequireJS, whose global define() causes Mermaid's UMD build
// to register as an AMD module instead of setting window.mermaid — so
// auto-init never fires. We work around this by injecting a type="module"
// script, which runs in the ES-module loader and is unaffected by RequireJS.
(function () {
  var s = document.createElement('script');
  s.type = 'module';
  s.textContent = [
    'import mermaid from "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs";',
    'mermaid.initialize({ startOnLoad: true, theme: "default" });',
  ].join('\n');
  document.head.appendChild(s);
})();
