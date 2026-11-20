"""NGL structure view with soft lighting and exact, pre-resolved atom colors."""

from __future__ import annotations

import json

from data_access import UNMAPPED_COLOR
from appearance import colors

NGL_URL = "https://cdn.jsdelivr.net/npm/ngl@2.4.0/dist/ngl.js"


def structure_payload(pdb_text: str, styles: list[dict]) -> dict:
    # Atom serials were resolved against the selected PDB model in mapped_colors.
    # Never substitute NGL atom indices or infer residue offsets here.
    colors = {}
    for style in styles:
        for serial in style["selection"]["serial"]:
            colors[str(int(serial))] = int(style["color"].lstrip("#"), 16)
    return {"pdb": pdb_text, "colors": colors, "unmapped": int(UNMAPPED_COLOR[1:], 16)}


def structure_html(pdb_text: str, styles: list[dict], *, dark: bool = False) -> str:
    # Escape HTML-significant characters even inside a JSON script element.
    payload = json.dumps(structure_payload(pdb_text, styles), separators=(",", ":"))
    payload = payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    palette = colors(dark)
    rendered = _TEMPLATE.replace("__NGL_URL__", NGL_URL).replace("__PAYLOAD__", payload)
    for key in ("viewer", "surface", "text", "border", "hover", "accent"):
        rendered = rendered.replace("__" + key.upper() + "__", palette[key])
    return rendered


_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<style>
html,body {margin:0;height:100%;font:13px system-ui,sans-serif;color:__TEXT__;background:__VIEWER__;}
.viewer {position:relative;height:100%;box-sizing:border-box;border:1px solid __BORDER__;border-radius:14px;overflow:hidden;}
#viewport {width:100%;height:100%;}
.toolbar {position:absolute;right:12px;top:12px;display:flex;gap:6px;z-index:2;}
button {background:__SURFACE__;border:1px solid __BORDER__;color:__TEXT__;padding:6px 10px;border-radius:7px;cursor:pointer;font:inherit;}
button:hover {background:__HOVER__;} button:focus-visible {outline:2px solid __ACCENT__;outline-offset:2px;}
#status {position:absolute;left:18px;bottom:14px;max-width:80%;line-height:1.5;pointer-events:none;}
</style></head><body>
<div class="viewer" aria-label="Interactive protein structure colored by protection factor">
  <div id="viewport"></div>
  <div class="toolbar"><button id="reset" disabled>Reset view</button></div>
  <div id="status" role="status">Loading structure…</div>
</div>
<script id="structure-data" type="application/json">__PAYLOAD__</script>
<script src="__NGL_URL__" onerror="document.getElementById('status').textContent='The structure viewer could not load. Check your connection and refresh.'"></script>
<script>
(async function () {
  const status = document.getElementById('status');
  if (!window.NGL) return;
  let stage;
  try {
    const data = JSON.parse(document.getElementById('structure-data').textContent);
    stage = new NGL.Stage('viewport', {
      backgroundColor: '__VIEWER__', cameraType: 'orthographic', quality: 'high',
      sampleLevel: 2, lightColor: '#ffffff', lightIntensity: 0.4,
      ambientColor: '#ffffff', ambientIntensity: 0.45,
      tooltip: true, workerDefault: false
    });
    if (!stage.viewer.renderer) throw new Error('WebGL is unavailable');
    const scheme = NGL.ColormakerRegistry.addScheme(function () {
      this.atomColor = function (atom) {
        return data.colors[String(atom.serial)] ?? data.unmapped;
      };
    });
    const component = await stage.loadFile(new Blob([data.pdb], {type: 'text/plain'}), {
      ext: 'pdb', defaultRepresentation: false, useWorker: false
    });
    component.addRepresentation('cartoon', {
      colorScheme: scheme, quality: 'high', smoothSheet: true,
      metalness: 0, roughness: 0.85, useInteriorColor: false,
      diffuseInterior: true, interiorDarkening: 0
    });
    component.addRepresentation('ball+stick', {
      sele: 'hetero and not (water or ion) and not polymer', colorScheme: 'element',
      radiusScale: 0.6, metalness: 0, roughness: 0.85
    });
    component.autoView(0);
    status.textContent = '';
    document.getElementById('reset').disabled = false;
    document.getElementById('reset').onclick = () => component.autoView(250);
    new ResizeObserver(() => stage.handleResize()).observe(document.getElementById('viewport'));
    window.addEventListener('pagehide', () => stage.dispose(), {once: true});
  } catch (error) {
    console.error('PFBoost structure viewer:', error);
    status.textContent = 'The structure could not be displayed. Check that WebGL is enabled, then refresh.';
    if (stage) stage.dispose();
  }
})();
</script></body></html>
"""
