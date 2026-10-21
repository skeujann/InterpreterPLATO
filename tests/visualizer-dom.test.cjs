'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { createRequire } = require('node:module');
const requireViewer = createRequire(path.join(__dirname, '../visualizer/package.json'));
const { JSDOM } = requireViewer('jsdom');
const fixture = require('../visualizer/example.js');
const root = path.join(__dirname, '../visualizer');
const tick = () => new Promise(resolve => setTimeout(resolve, 0));
const copy = () => JSON.parse(JSON.stringify(fixture));

async function harness() {
  const dom = new JSDOM(fs.readFileSync(path.join(root, 'index.html'), 'utf8'), { runScripts: 'outside-only', url: 'file:///InterpreterPLATO/visualizer/index.html' });
  const { window } = dom;
  const errors = [];
  window.structuredClone = globalThis.structuredClone;
  window.addEventListener('error', event => { errors.push(event.error || event.message); event.preventDefault(); });
  if (window.document.readyState === 'loading') await new Promise(resolve => window.document.addEventListener('DOMContentLoaded', resolve, { once: true }));
  const instances = [];
  window.eval(fs.readFileSync(path.join(root, 'vendor/cytoscape.min.js'), 'utf8'));
  const browserCytoscape = window.cytoscape;
  window.cytoscape = options => {
    // Load both the adapter and Cytoscape in the same DOM realm, as browsers do.
    const configuration = new window.Object();
    Object.assign(configuration, options, { container: undefined, headless: true, styleEnabled: false });
    const cy = browserCytoscape(configuration);
    instances.push(cy); return cy;
  };
  window.fetch = () => { throw new Error('Network access is forbidden in viewer tests'); };
  window.alert = () => { throw new Error('Untrusted text was executed'); };
  for (const source of ['example.js', 'graph-model.js', 'app.js']) window.eval(fs.readFileSync(path.join(root, source), 'utf8'));
  await tick();
  return { window, document: window.document, cy: () => instances[instances.length - 1],
    close: () => { instances.forEach(cy => cy.destroy()); window.close(); assert.deepEqual(errors, [], 'No unhandled UI event errors'); } };
}
async function upload(h, payload, name = 'explanation.json') {
  const text = typeof payload === 'string' ? payload : JSON.stringify(payload);
  const input = h.document.getElementById('file-input');
  Object.defineProperty(input, 'files', { configurable: true, value: [{ name, size: text.length, text: async () => text }] });
  input.dispatchEvent(new h.window.Event('change', { bubbles: true })); await tick();
}

test('UI starts with an explicitly non-computational synthetic example and hides report structure', async () => {
  const h = await harness();
  try {
    assert.match(h.document.getElementById('mode-title').textContent, /非真实计算结果/);
    assert.match(h.document.getElementById('case-title').textContent, /Demo-A \+ Demo-B/);
    assert.match(h.document.getElementById('graph-counts').textContent, /6 条 PPI/);
    assert.equal(h.cy().nodes().filter(node => node.data('group') === 'report' && !node.hasClass('filtered')).length, 0);
    assert.equal(h.document.getElementById('error-banner').hidden, true);
  } finally { h.close(); }
});
test('path buttons, search, filters, inspector and layouts act on the loaded graph', async () => {
  const h = await harness();
  try {
    h.document.querySelector('.path-button').click();
    assert.equal(h.cy().edges('.highlighted').length, 3);
    assert.equal(h.cy().edges('.highlighted').every(edge => edge.data('type') === 'PPI'), true);
    const search = h.document.getElementById('node-search'); search.value = 'SYN_UNMAPPED'; search.dispatchEvent(new h.window.Event('input'));
    assert.equal(h.cy().nodes('.highlighted').length, 1);
    search.dispatchEvent(new h.window.KeyboardEvent('keydown', { key: 'Enter' }));
    assert.equal(h.document.getElementById('detail-title').textContent, 'SYN_UNMAPPED');
    const filter = h.document.querySelector('[data-node-group="gene"]'); filter.checked = false; filter.dispatchEvent(new h.window.Event('change'));
    assert.equal(h.cy().edges().filter(edge => edge.data('type') === 'PPI' && !edge.hasClass('filtered')).length, 0);
    filter.checked = true; filter.dispatchEvent(new h.window.Event('change'));
    h.cy().nodes()[0].emit('tap'); assert.equal(h.document.getElementById('detail-json').hidden, false);
    const layout = h.document.getElementById('layout-select');
    for (const name of ['breadthfirst', 'concentric', 'cose', 'grid']) { layout.value = name; layout.dispatchEvent(new h.window.Event('change')); }
    h.document.getElementById('clear-selection').click(); assert.equal(h.cy().elements('.highlighted').length, 0);
  } finally { h.close(); }
});
test('uploads stay local, malformed files do not replace the current record, and labels are never injected as HTML', async () => {
  const h = await harness();
  try {
    const payload = copy(); payload.drug1_name = '<img src=x onerror=alert(1)>';
    payload.evidence_graphs.explanation_evidence_graph.nodes[0].name = payload.drug1_name;
    await upload(h, payload, '<b>local.json</b>');
    assert.match(h.document.getElementById('case-title').textContent, /<img/);
    assert.equal(h.document.getElementById('source-name').textContent, '<b>local.json</b>');
    assert.equal(h.document.querySelector('img'), null); assert.equal(h.document.querySelector('#source-name b'), null);
    h.cy().nodes()[0].emit('tap'); assert.match(h.document.getElementById('detail-json').textContent, /<img/);
    await upload(h, '{invalid'); assert.equal(h.document.getElementById('error-banner').hidden, false);
    assert.match(h.document.getElementById('case-title').textContent, /<img/);
  } finally { h.close(); }
});
test('batch selection retains failed records and displays their diagnostics without fabricated nodes', async () => {
  const h = await harness();
  try {
    const invalid = copy(); invalid.drug1_name = 'Invalid graph'; invalid.evidence_graphs.explanation_evidence_graph.edges[0].target = 'missing';
    await upload(h, [copy(), { drug1: 'Failed A', drug2: 'Failed B', error: 'No eligible human targets' }, invalid]);
    const select = h.document.getElementById('result-select'); assert.equal(select.hidden, false);
    select.value = '1'; select.dispatchEvent(new h.window.Event('change'));
    assert.equal(h.cy().nodes().length, 0); assert.equal(h.document.getElementById('graph-empty').hidden, false);
    assert.equal(h.document.getElementById('export-png').disabled, true);
    assert.match(h.document.getElementById('diagnostic-summary').textContent, /No eligible human targets/);
    assert.equal(h.document.getElementById('diagnostics-panel').hidden, false);
    select.value = '2'; select.dispatchEvent(new h.window.Event('change'));
    assert.equal(h.cy().nodes().length, 0); assert.match(h.document.getElementById('case-title').textContent, /Invalid graph/);
    assert.match(h.document.getElementById('diagnostic-summary').textContent, /缺失节点/);
    assert.equal(h.document.getElementById('detail-json').hidden, true);
  } finally { h.close(); }
});
test('real-input test mode is not relabeled synthetic and does not acquire a validation claim', async () => {
  const h = await harness();
  try {
    const payload = copy(); payload.generated_from_test_fixture = false; payload.metadata = {}; payload.run_mode = 'test';
    await upload(h, payload);
    assert.match(h.document.getElementById('mode-title').textContent, /测试／演示模式/);
    assert.doesNotMatch(h.document.getElementById('mode-title').textContent, /合成/);
    assert.match(h.document.getElementById('mode-description').textContent, /不独立核验/);
  } finally { h.close(); }
});
test('workspace panels preserve selections and filters and open details when requested', async () => {
  const h = await harness();
  try {
    const get = id => h.document.getElementById(id);
    assert.equal(get('inspector-panel').hidden, true);
    assert.equal(get('path-count').textContent, String(h.document.querySelectorAll('.path-button').length));
    get('tab-filters').click();
    assert.equal(get('paths-panel').hidden, true); assert.equal(get('filters-panel').hidden, false);
    const filter = h.document.querySelector('[data-node-group="gene"]'); filter.checked = false; filter.dispatchEvent(new h.window.Event('change'));
    get('toggle-sidebar').click();
    assert.equal(get('explore-panel').hidden, true); assert.equal(get('toggle-sidebar').getAttribute('aria-expanded'), 'false');
    get('toggle-sidebar').click();
    assert.equal(filter.checked, false); assert.equal(get('filters-panel').hidden, false);
    filter.checked = true; filter.dispatchEvent(new h.window.Event('change'));
    get('tab-paths').click(); h.document.querySelector('.path-button').click();
    assert.equal(get('inspector-panel').hidden, false); assert.equal(get('toggle-inspector').getAttribute('aria-expanded'), 'true');
    assert.equal(h.document.querySelector('.path-button').getAttribute('aria-pressed'), 'true');
    const selected = h.cy().elements('.highlighted').map(item => item.id());
    get('close-inspector').click();
    assert.equal(get('inspector-panel').hidden, true); assert.equal(h.document.activeElement, get('toggle-inspector'));
    assert.deepEqual(h.cy().elements('.highlighted').map(item => item.id()), selected);
    get('toggle-inspector').click(); get('inspector-panel').dispatchEvent(new h.window.KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    assert.equal(get('inspector-panel').hidden, true);
    h.cy().nodes()[0].emit('tap'); assert.equal(get('inspector-panel').hidden, false);
    await upload(h, copy()); assert.equal(get('inspector-panel').hidden, true);
  } finally { h.close(); }
});

test('tabs support arrow-key navigation and zoom controls respect viewport limits', async () => {
  const h = await harness();
  try {
    const get = id => h.document.getElementById(id);
    get('tab-paths').focus(); get('tab-paths').dispatchEvent(new h.window.KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    assert.equal(h.document.activeElement, get('tab-filters')); assert.equal(get('tab-filters').tabIndex, 0); assert.equal(get('tab-paths').tabIndex, -1);
    get('tab-filters').dispatchEvent(new h.window.KeyboardEvent('keydown', { key: 'Home', bubbles: true }));
    assert.equal(get('paths-panel').hidden, false);
    get('toggle-inspector').click(); get('tab-details').dispatchEvent(new h.window.KeyboardEvent('keydown', { key: 'End', bubbles: true }));
    assert.equal(get('diagnostics-panel').hidden, false); assert.equal(h.document.activeElement, get('tab-diagnostics'));
    h.cy().zoom(1); get('zoom-in').click(); assert.equal(h.cy().zoom(), 1.25);
    get('zoom-out').click(); assert.equal(h.cy().zoom(), 1);
    for (let i = 0; i < 30; i += 1) get('zoom-in').click();
    assert.equal(h.cy().zoom(), h.cy().maxZoom());
    for (let i = 0; i < 60; i += 1) get('zoom-out').click();
    assert.equal(h.cy().zoom(), h.cy().minZoom());
    await upload(h, { drug1: 'A', drug2: 'B', error: 'No evidence graph' });
    for (const id of ['zoom-in', 'zoom-out', 'fit-graph', 'export-png']) assert.equal(get(id).disabled, true);
  } finally { h.close(); }
});

test('PNG export includes the synthetic/non-computational watermark before download (canvas stub)', async () => {
  const h = await harness();
  try {
    const drawn = [], downloads = [];
    h.cy().png = () => 'data:image/png;base64,canvas-test-only';
    h.window.Image = class {
      constructor() { this.width = 500; this.height = 300; }
      set src(value) { queueMicrotask(() => this.onload()); }
    };
    h.window.HTMLCanvasElement.prototype.getContext = () => ({ fillRect() {}, drawImage() {}, fillText(value) { drawn.push(value); } });
    h.window.HTMLCanvasElement.prototype.toBlob = function (callback) { callback(new h.window.Blob(['canvas-test-only'])); };
    h.window.URL.createObjectURL = () => 'blob:local-test'; h.window.URL.revokeObjectURL = () => {};
    h.window.HTMLAnchorElement.prototype.click = function () { downloads.push(this.download); };
    h.document.getElementById('export-png').click(); await tick();
    assert.match(drawn.join(' '), /非真实计算结果/); assert.match(drawn.join(' '), /未独立核验/);
    assert.deepEqual(downloads, ['InterpreterPLATO-SYNTHETIC-network.png']);
  } finally { h.close(); }
});
