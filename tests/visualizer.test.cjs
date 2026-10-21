'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const vm = require('node:vm');
const cytoscape = require('../visualizer/vendor/cytoscape.min.js');
const G = require('../visualizer/graph-model.js');
const fixture = require('../visualizer/example.js');
const { build } = require('../visualizer/build.cjs');
const copy = () => JSON.parse(JSON.stringify(fixture));

test('preserves exactly the supplied graph without changing the input', () => {
  const payload = copy(), before = JSON.stringify(payload), model = G.build(payload);
  const input = payload.evidence_graphs.explanation_evidence_graph;
  assert.equal(model.elements.filter(item => item.group === 'nodes').length, input.nodes.length);
  assert.equal(model.elements.filter(item => item.group === 'edges').length, input.edges.length);
  assert.equal(model.ppiCount, 6);
  assert.equal(JSON.stringify(payload), before);
  assert.equal(new Set(model.elements.map(item => item.data.id)).size, model.elements.length);
});
test('Cytoscape loads the adapter output and path highlights contain actual PPI only', () => {
  const model = G.build(copy());
  const cy = cytoscape({ headless: true, elements: model.elements });
  assert.equal(cy.nodes().length, 14); assert.equal(cy.edges().length, 19);
  for (const row of model.paths) {
    assert.equal(row.valid, true);
    row.edgeIds.forEach(id => { const edge = cy.getElementById(id); assert.equal(edge.data('type'), 'PPI'); });
    row.nodeIds.forEach(id => assert.equal(cy.getElementById(id).isNode(), true));
  }
  cy.destroy();
});
test('never reconstructs a missing graph from targets or path lists', () => {
  const payload = copy(); delete payload.evidence_graphs;
  const result = G.build(payload);
  assert.equal(result.elements.length, 0); assert.equal(result.paths.length, 0);
  assert.match(result.warnings[0], /不会/);
});
test('a path missing a recorded PPI is disabled instead of acquiring a fabricated edge', () => {
  const payload = copy(), graph = payload.evidence_graphs.explanation_evidence_graph;
  graph.edges = graph.edges.filter(edge => !(edge.type === 'PPI' && edge.source === 'gene:SYN_A'));
  const model = G.build(payload), cross = model.paths.find(row => row.kind === 'cross' && !row.zeroHop);
  assert.equal(cross.valid, false); assert.match(cross.missing.join(' '), /缺少 PPI/);
  assert.equal(model.elements.filter(row => row.group === 'edges').length, graph.edges.length);
});
test('zero-hop membership highlights one node without inventing a self edge', () => {
  const model = G.build(copy()), zero = model.paths.find(row => row.zeroHop);
  assert.equal(zero.valid, true); assert.equal(zero.nodeIds.length, 1); assert.equal(zero.edgeIds.length, 0);
  assert.equal(model.elements.filter(row => row.group === 'edges' && row.data.source === row.data.target).length, 0);
});
test('all visual relationship categories suppress every arrow', () => {
  const base = G.styles.find(row => row.selector === 'edge').style;
  for (const property of ['source-arrow-shape', 'target-arrow-shape', 'mid-source-arrow-shape', 'mid-target-arrow-shape']) assert.equal(base[property], 'none');
  for (const row of G.styles) for (const [property, value] of Object.entries(row.style)) if (property.endsWith('arrow-shape')) assert.equal(value, 'none');
  assert.notEqual(G.COLORS.ppi, G.COLORS.disease); assert.notEqual(G.COLORS.ppi, G.COLORS.report);
  assert.equal(G.relationship('belongs_to_disease_module'), 'disease');
  assert.equal(G.relationship('summarizes_combination_member'), 'report');
  assert.equal(G.relationship('arbitrary_biological_arrow'), 'other');
});
test('rejects dangling edges, duplicate nodes, directed PPI and non-gene PPI', () => {
  let payload = copy(); payload.evidence_graphs.explanation_evidence_graph.edges[0].target = 'missing'; assert.throws(() => G.build(payload), /缺失节点/);
  payload = copy(); payload.evidence_graphs.explanation_evidence_graph.nodes.push(payload.evidence_graphs.explanation_evidence_graph.nodes[0]); assert.throws(() => G.build(payload), /重复/);
  payload = copy(); payload.evidence_graphs.explanation_evidence_graph.edges[0].directed = true; assert.throws(() => G.build(payload), /无向/);
  payload = copy(); payload.evidence_graphs.explanation_evidence_graph.edges[0].source = 'drug:Demo-A'; assert.throws(() => G.build(payload), /基因/);
});
test('distinguishes synthetic fixtures from real-input test mode without granting scientific validation', () => {
  assert.equal(G.isSynthetic(fixture), true); assert.match(G.modeLabel(fixture), /非真实计算结果/);
  assert.equal(G.isSynthetic({ run_mode: 'test', generated_from_test_fixture: false }), false);
  assert.match(G.modeLabel({ run_mode: 'test' }), /测试/);
  assert.match(G.modeLabel({ run_mode: 'research_hypothesis', formal_result_eligible: true }), /非因果/);
  assert.match(G.modeLabel({}), /未知/);
  assert.equal(fixture.metadata.calculation_performed, false);
  assert.equal(fixture.network_statistics.valid, false);
});
test('preserves untrusted labels as data and supports IDs that resemble object keys', () => {
  const payload = copy(), graph = payload.evidence_graphs.explanation_evidence_graph;
  graph.nodes[0].name = '</script><img src=x onerror=alert(1)>';
  graph.nodes.push({ id: '__proto__', type: 'gene', name: 'constructor' });
  const model = G.build(payload);
  assert.equal(model.elements[0].data.label, graph.nodes[0].name);
  assert.equal(model.elements.filter(row => row.data.originalId === '__proto__').length, 1);
});
test('supports explanation arrays and retains failed records, but rejects cache/manifests and oversized inputs', () => {
  const failed = { drug1: 'A', drug2: 'B', error: 'No eligible targets' };
  assert.equal(G.records([fixture, failed]).length, 2); assert.equal(G.records({ results: [fixture] }).length, 1);
  assert.match(G.build(failed).warnings[0], /未生成解释/);
  assert.throws(() => G.records({ cache_schema_version: 3 }), /解释结果 JSON/);
  assert.throws(() => G.records([])); assert.throws(() => G.records(Array(101).fill(fixture)));
  const payload = copy(); payload.evidence_graphs.explanation_evidence_graph.nodes = Array(G.LIMITS.nodes + 1).fill({ id: 'x' });
  assert.throws(() => G.build(payload), /上限/);
});
test('unknown edge types remain unclassified and never become PPI', () => {
  const payload = copy(); payload.evidence_graphs.explanation_evidence_graph.edges[0].type = 'custom_relation';
  const model = G.build(payload); assert.equal(model.ppiCount, 5); assert.match(model.warnings[0], /未分类/);
});
test('the optional single-file distribution embeds every asset, the license and syntactically valid scripts', () => {
  const temporary = fs.mkdtempSync(path.join(os.tmpdir(), 'plato-visualizer-'));
  try {
    const output = build(path.join(temporary, 'viewer.html')), html = fs.readFileSync(output, 'utf8');
    assert.doesNotMatch(html, /<script[^>]+src=/i); assert.doesNotMatch(html, /<link[^>]+stylesheet/i);
    assert.match(html, /Permission is hereby granted/); assert.match(html, /Cytoscape\.js 3\.34\.3/);
    const scripts = [...html.matchAll(/<script>\n([\s\S]*?)\n<\/script>/g)];
    assert.equal(scripts.length, 5); scripts.forEach(([, source]) => new vm.Script(source));
  } finally {
    assert.equal(path.dirname(path.resolve(temporary)), path.resolve(os.tmpdir()));
    assert.ok(path.basename(temporary).startsWith('plato-visualizer-'));
    fs.rmSync(temporary, { recursive: true, force: true });
  }
});
