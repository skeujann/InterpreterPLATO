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

async function harness(example = 'synthetic', runBundle, styleEnabled = false) {
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
    Object.assign(configuration, options, { container: undefined, headless: true, styleEnabled });
    const cy = browserCytoscape(configuration);
    instances.push(cy); return cy;
  };
  window.fetch = () => { throw new Error('Network access is forbidden in viewer tests'); };
  window.alert = () => { throw new Error('Untrusted text was executed'); };
  if (runBundle) {
    const value = window.JSON.parse(JSON.stringify(runBundle));
    if (Array.isArray(value)) { window.InterpreterRunExamples = value; window.InterpreterRunExample = value[0]; }
    else window.InterpreterRunExample = value;
  }
  for (const source of ['example.js', 'reference.js', 'graph-model.js', 'app.js']) window.eval(fs.readFileSync(path.join(root, source), 'utf8'));
  await tick();
  if (example === 'synthetic') {
    const select = window.document.getElementById('example-select'); select.value = 'synthetic'; select.dispatchEvent(new window.Event('change')); await tick();
  }
  return { window, document: window.document, cy: () => instances[instances.length - 1],
    close: () => { instances.forEach(cy => cy.destroy()); window.close(); assert.deepEqual(errors, [], 'No unhandled UI event errors'); } };
}
async function upload(h, payload, name = 'explanation.json') {
  const text = typeof payload === 'string' ? payload : JSON.stringify(payload);
  const input = h.document.getElementById('file-input');
  Object.defineProperty(input, 'files', { configurable: true, value: [{ name, size: text.length, text: async () => text }] });
  input.dispatchEvent(new h.window.Event('change', { bubbles: true })); await tick();
}

test('synthetic example remains explicitly non-computational and hides report structure', async () => {
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
    h.cy().nodes()[0].emit('mouseover');
    assert.equal(h.document.querySelector('#graph-tooltip img'), null); assert.match(h.document.getElementById('graph-tooltip').textContent, /<img/);
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

test('default reference case has a guided layout and highlights only supplied relationship types', async () => {
  const h = await harness('reference');
  try {
    const get = id => h.document.getElementById(id);
    assert.equal(get('case-guide').hidden, false); assert.equal(get('inspector-panel').hidden, true);
    assert.match(get('case-title').textContent, /Dabrafenib \+ Trametinib/);
    assert.match(get('mode-title').textContent, /文献导览.*非模型计算/);
    assert.equal(h.cy().nodes().length, 11); assert.equal(h.cy().edges().length, 20);
    const shown = item => !item.hasClass('filtered') && !item.hasClass('focus-hidden');
    assert.equal(h.cy().nodes().filter(shown).length, 6);
    assert.equal(h.cy().edges().filter(shown).length, 6);
    assert.match(get('graph-counts').textContent, /2 条 PPI/); assert.equal(get('layout-select').value, 'preset');
    assert.equal(h.document.querySelectorAll('.path-button').length, 3);
    assert.equal(get('view-focus').getAttribute('aria-pressed'), 'true');
    assert.equal(get('focus-stage-labels').hidden, false);
    const endpoint = h.cy().nodes('.reference-context');
    assert.equal(endpoint.length, 1); assert.equal(endpoint.data('originalId'), 'disease_module:literature-context');
    assert.match(endpoint.data('displayLabel'), /疾病模块.*黑色素瘤/);
    assert.match(endpoint.data('displayLabel'), /文献背景/);
    const anchor = h.cy().nodes('.reference-anchor');
    assert.equal(anchor.data('originalId'), 'gene:BRAF');
    assert.ok(h.cy().nodes().filter(node => shown(node) && node.id() !== endpoint.id()).every(node => node.position('x') < endpoint.position('x')));
    const ids = h.cy().elements().map(item => item.id());
    get('view-full').click();
    assert.equal(h.cy().nodes().filter(shown).length, 11); assert.equal(h.cy().edges().filter(shown).length, 20);
    assert.match(get('graph-counts').textContent, /9 条 PPI/); assert.equal(h.document.querySelectorAll('.path-button').length, 6);
    assert.equal(get('focus-stage-labels').hidden, true);
    get('view-focus').click();
    assert.equal(h.cy().nodes().filter(shown).length, 6); assert.deepEqual(h.cy().elements().map(item => item.id()), ids);
    const paths = [...h.document.querySelectorAll('.path-button')];
    paths[1].click();
    assert.equal(h.cy().edges('.highlighted').length, 3);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('type') === 'PPI').length, 1);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('type') === 'known_target').length, 1);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('category') === 'disease').length, 1);
    assert.match(get('detail-title').textContent, /Trametinib.*疾病模块/);
    paths[0].click(); assert.equal(h.cy().edges('.highlighted').length, 2);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('type') === 'known_target').length, 1);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('category') === 'disease').length, 1);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('type') === 'PPI').length, 0);
    assert.equal(h.document.querySelectorAll('.path-button:disabled').length, 0);
    assert.equal(get('case-guide').open, false);
    get('case-guide').querySelector('summary').click(); assert.equal(get('case-guide').open, true);
    h.document.querySelector('[data-tour="targets"]').click();
    assert.equal(h.cy().edges('.highlighted').length, 3);
    assert.ok(h.cy().edges('.highlighted').every(edge => edge.data('type') === 'known_target'));
    h.document.querySelector('[data-tour="ppi"]').click();
    assert.equal(h.cy().edges('.highlighted').length, 2);
    assert.ok(h.cy().edges('.highlighted').every(edge => edge.data('type') === 'PPI'));
    assert.match(get('tour-status').textContent, /不是因果/);
    h.document.querySelector('[data-tour="context"]').click();
    assert.equal(h.cy().edges('.highlighted').length, 8);
    assert.ok(h.cy().edges('.highlighted').every(edge => ['disease', 'pathway'].includes(edge.data('category'))));
    assert.equal(h.cy().nodes().filter(node => !node.hasClass('filtered')).length, 11);
    assert.match(get('detail-description').textContent, /单基因/);
    get('show-reference-sources').click();
    assert.equal(get('reference-sources').hidden, false); assert.match(get('detail-title').textContent, /来源/);
    for (const link of h.document.querySelectorAll('#reference-sources a')) { assert.match(link.href, /^https:\/\//); assert.equal(link.rel, 'noopener noreferrer'); }
    h.document.querySelector('[data-tour="overview"]').click();
    assert.equal(h.cy().elements('.highlighted').length, 0); assert.equal(get('inspector-panel').hidden, true);
    assert.equal(h.cy().nodes().length, 11); assert.equal(h.cy().edges().length, 20);
    Object.defineProperty(get('cy'), 'clientWidth', { configurable: true, value: 390 });
    h.window.dispatchEvent(new h.window.Event('resize'));
    assert.equal(get('focus-stage-labels').hidden, true);
    assert.ok(h.cy().nodes('[type="drug"]').every(node => node.position('y') < endpoint.position('y')));
    assert.equal(h.cy().nodes('.compact-reference').length, 11);
    assert.deepEqual(h.cy().elements().map(item => item.id()), ids, 'Narrow-screen layout preserves the complete graph');
  } finally { h.close(); }
});

test('uploaded files cannot activate curated layout, source links or reference tour by spoofing metadata', async () => {
  const h = await harness('reference');
  try {
    const payload = copy(); payload.metadata.reference_example = true;
    payload.metadata.visualization_positions = { 'drug:Demo-A': [1000000, -1000000] };
    payload.metadata.references = [{ title: 'untrusted', url: 'javascript:alert(1)' }];
    await upload(h, payload);
    assert.equal(h.document.getElementById('case-guide').hidden, true);
    assert.equal(h.document.getElementById('reference-sources').hidden, true);
    assert.equal(h.document.getElementById('reference-layout').disabled, true);
    assert.equal(h.document.getElementById('layout-select').value, 'cose');
    assert.equal(h.document.getElementById('reference-view-switch').hidden, true);
    assert.equal(h.document.getElementById('focus-stage-labels').hidden, true);
    assert.equal(h.cy().nodes('.reference-anchor, .reference-context, .reference-drug, .reference-target').length, 0);
    assert.match(h.document.getElementById('hero-endpoint').textContent, /JSON/);
    assert.equal(h.document.querySelector('a[href^="javascript:"]'), null);
    const model = h.cy().elements().map(item => item.id());
    h.document.querySelector('[data-tour="ppi"]').click();
    assert.deepEqual(h.cy().elements().map(item => item.id()), model);
    h.document.getElementById('example-select').value = 'reference'; h.document.getElementById('load-example').click();
    assert.equal(h.document.getElementById('case-guide').hidden, false); assert.equal(h.cy().nodes().length, 11);
  } finally { h.close(); }
});

test('frozen reference JSON reloads without changing its reference-only semantics', async () => {
  const h = await harness('reference');
  try {
    const contents = fs.readFileSync(path.join(root, 'examples/dabrafenib-trametinib.json'), 'utf8');
    const payload = JSON.parse(contents);
    assert.equal(payload.run_mode, 'reference_example'); assert.equal(payload.predicted_synergy, null);
    await upload(h, contents, 'dabrafenib-trametinib.reference.json');
    assert.equal(h.cy().nodes().length, 11); assert.equal(h.cy().edges().length, 20);
    assert.equal(h.document.querySelectorAll('.path-button:disabled').length, 0);
    assert.match(h.document.getElementById('path-list').textContent, /文献背景/);
    assert.match(h.document.getElementById('mode-title').textContent, /非模型计算结果/);
    assert.equal(h.document.getElementById('case-guide').hidden, true);
    assert.equal(h.document.getElementById('reference-layout').disabled, true);
  } finally { h.close(); }
});

test('embedded run preserves actual statistics, complete graph, path types and upload isolation', async () => {
  const { withRunFile } = require('./run-example-fixture.cjs');
  const { prepareRun } = require('../visualizer/prepare-run.cjs');
  const bundle = withRunFile(({ file }) => prepareRun(file));
  const h = await harness('run', bundle, true);
  try {
    const get = id => h.document.getElementById(id), shown = item => !item.hasClass('filtered') && !item.hasClass('focus-hidden');
    assert.equal(get('example-select').value, 'run'); assert.equal(get('run-summary').hidden, false);
    assert.equal(get('mode-banner').hidden, true);
    assert.equal(get('mode-title').textContent, ''); assert.equal(get('mode-description').textContent, '');
    assert.match(get('run-verdict').textContent, /未获显著支持/);
    assert.match(get('run-metrics').textContent, /0\.8000/); assert.match(get('run-metrics').textContent, /1\.0000/);
    assert.match(get('run-metrics').textContent, /调用方提供/); assert.match(get('run-metrics').textContent, /12\.3457/);
    assert.equal(h.cy().nodes().filter(shown).length, 6); assert.equal(h.cy().edges().filter(shown).length, 4);
    assert.ok(h.cy().nodes().filter(shown).every(node => node.style('label') !== ''));
    assert.equal(get('case-guide').hidden, true); assert.equal(get('run-guide').hidden, false);
    assert.equal(h.document.querySelectorAll('.path-button').length, 2);
    const terminal = h.cy().nodes('.run-context');
    assert.equal(terminal.length, 1); assert.equal(terminal.hasClass('focus-hidden'), true);
    assert.doesNotMatch(get('path-list').textContent, /fixture_context/);
    assert.equal(get('focus-stage-labels').children.length, 3);
    h.document.querySelector('.path-button').click();
    assert.deepEqual([...h.cy().edges('.highlighted')].map(edge => edge.data('category')).sort(), ['ppi', 'target']);
    const ids = h.cy().elements().map(item => item.id());
    get('view-full').click(); assert.equal(h.cy().nodes().filter(shown).length, 9); assert.equal(h.cy().edges().filter(shown).length, 8);
    assert.equal(shown(terminal), true);
    assert.ok(h.cy().nodes().every(node => node.style('label') === ''));
    h.document.querySelector('.path-button').click();
    const selected = h.cy().elements('.highlighted').map(item => item.id());
    const hovered = h.cy().nodes().filter(node => node.data('originalId') === 'gene:A0');
    hovered.emit('mouseover');
    assert.equal(get('graph-tooltip').hidden, false); assert.match(get('graph-tooltip').textContent, /A0/);
    assert.match(get('graph-tooltip').textContent, /3 个可见相邻节点 · 3 条连接/);
    assert.deepEqual(h.cy().edges('.hover-connection').map(edge => edge.id()).sort(), hovered.connectedEdges().filter(shown).map(edge => edge.id()).sort());
    hovered.emit('mouseout'); assert.equal(get('graph-tooltip').hidden, true);
    assert.equal(h.cy().elements('.hover-neighbor, .hover-connection, .hovered').length, 0);
    assert.deepEqual(h.cy().elements('.highlighted').map(item => item.id()), selected);
    const memberships = h.document.querySelector('[data-edge-group="membership"]'); memberships.checked = false; memberships.dispatchEvent(new h.window.Event('change'));
    hovered.emit('mouseover'); assert.match(get('graph-tooltip').textContent, /2 个可见相邻节点 · 2 条连接/);
    assert.ok(h.cy().edges('.hover-connection').every(shown));
    get('view-focus').click(); assert.deepEqual(h.cy().elements().map(item => item.id()), ids);
    assert.equal(get('graph-tooltip').hidden, true); assert.ok(h.cy().nodes().filter(shown).every(node => node.style('label') !== ''));
    get('show-run-details').click(); assert.match(get('detail-json').textContent, /false/); assert.match(get('detail-json').textContent, new RegExp(bundle.source.sha256));
    const payload = JSON.parse(JSON.stringify(bundle.payload)); payload.metadata = { ...bundle.presentation };
    await upload(h, payload);
    assert.equal(get('run-summary').hidden, true); assert.equal(get('reference-view-switch').hidden, true);
    assert.equal(get('mode-banner').hidden, false);
    assert.equal(h.cy().nodes('.run-drug, .run-context').length, 0);
    get('example-select').value = 'run'; get('load-example').click(); assert.equal(get('run-summary').hidden, false);
    get('example-select').value = 'reference'; get('load-example').click();
    assert.equal(get('run-summary').hidden, true); assert.match(get('focus-stage-labels').textContent, /BRAF/);
  } finally { h.close(); }
});

test('run suite switches case payloads and never promotes invalid statistics or zero-hop membership', async () => {
  const { withRunFile } = require('./run-example-fixture.cjs'), { prepareRun } = require('../visualizer/prepare-run.cjs');
  const bundles = withRunFile(({ file, payload }) => {
    const first = prepareRun(file);
    payload.audit.engineering_case.case_id = 'TEST-INVALID'; payload.cell_line = 'INVALID-CELL';
    Object.assign(payload.network_statistics, { valid: false, statistically_supported: true, quality_status: 'insufficient', p_A_adjusted: 0.00001 });
    payload.local_evidence.top_paths[0].path = ['A0']; payload.local_evidence.top_paths[0].length = 0;
    payload.evidence_graphs.explanation_evidence_graph.edges.push({ source: 'gene:A0', target: 'disease_module:fixture', type: 'belongs_to_disease_module' });
    fs.writeFileSync(file, JSON.stringify(payload)); const second = prepareRun(file);
    payload.audit.engineering_case.case_id = 'TEST-SINGLE-SIDE';
    Object.assign(payload.network_statistics, { valid: true, quality_status: 'limited', p_A_adjusted: 0.007, near_A: true, near_B: false });
    fs.writeFileSync(file, JSON.stringify(payload)); return [first, second, prepareRun(file)];
  });
  const h = await harness('run', bundles);
  try {
    const get = id => h.document.getElementById(id);
    const selectRun = value => { get('example-select').value = value; get('example-select').dispatchEvent(new h.window.Event('change')); };
    assert.deepEqual([...get('example-select').options].filter(option => option.value.startsWith('run')).map(option => option.value), ['run', 'run:1', 'run:2']);
    selectRun('run:1');
    assert.equal(get('example-select').value, 'run:1'); assert.match(get('case-meta').textContent, /INVALID-CELL/);
    assert.match(get('run-verdict').textContent, /不可评估/); assert.doesNotMatch(get('run-metrics').textContent, /0\.0000/);
    assert.match(get('run-metrics').textContent, /不可评估/);
    h.document.querySelector('.path-button').click();
    assert.equal(h.cy().edges('.highlighted').length, 1);
    assert.equal(h.cy().edges('.highlighted').filter(edge => edge.data('category') === 'ppi').length, 0);
    assert.match(get('search-status').textContent, /零跳/);
    selectRun('run:2');
    assert.match(get('run-verdict').textContent, /仅一侧/); assert.match(get('run-rule').textContent, /有限/);
    assert.match(get('run-metrics').textContent, /0\.0070/);
    selectRun('run');
    assert.match(get('run-verdict').textContent, /未获显著支持/); assert.match(get('run-metrics').textContent, /0\.8000/);
    assert.equal(get('example-select').value, 'run');
    const drawn = []; let pendingImage;
    h.cy().png = () => 'data:image/png;base64,test';
    h.window.Image = class { constructor() { pendingImage = this; this.width = 500; this.height = 300; } set src(value) {} };
    h.window.HTMLCanvasElement.prototype.getContext = () => ({ fillRect() {}, drawImage() {}, fillText(value) { drawn.push(value); } });
    h.window.HTMLCanvasElement.prototype.toBlob = callback => callback(new h.window.Blob(['test']));
    h.window.URL.createObjectURL = () => 'blob:local-export-test'; h.window.URL.revokeObjectURL = () => {};
    h.window.HTMLAnchorElement.prototype.click = () => {};
    get('export-png').click(); selectRun('run:1'); pendingImage.onload();
    assert.match(drawn.join(' '), /实际运行 TEST-CASE/);
    assert.doesNotMatch(drawn.join(' '), /TEST-INVALID/);
  } finally { h.close(); }
});

for (const example of ['synthetic', 'reference']) test('PNG export retains non-computational identity: ' + example + ' (canvas stub)', async () => {
  const h = await harness(example);
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
    assert.match(drawn.join(' '), example === 'reference' ? /文献导览.*非模型计算结果/ : /非真实计算结果/); assert.match(drawn.join(' '), /未独立核验/);
    assert.deepEqual(downloads, ['InterpreterPLATO-' + (example === 'reference' ? 'REFERENCE' : 'SYNTHETIC') + '-network.png']);
  } finally { h.close(); }
});
