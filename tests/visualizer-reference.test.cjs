'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const reference = require('../visualizer/reference.js');
const G = require('../visualizer/graph-model.js');
const pair = (a, b) => [a, b].sort().join('|');

test('reference module matches the frozen downloadable JSON and its byte hash', () => {
  const file = path.join(__dirname, '../visualizer/examples/dabrafenib-trametinib.json');
  const bytes = fs.readFileSync(file);
  assert.equal(bytes.includes(13), false, 'The pinned reference JSON must use LF on every platform');
  assert.equal(crypto.createHash('sha256').update(bytes).digest('hex'), fs.readFileSync(file + '.sha256', 'utf8').trim().split(/\s+/)[0]);
  assert.deepEqual(reference, JSON.parse(bytes.toString('utf8')));
  assert.equal(reference.run_mode, 'reference_example'); assert.equal(reference.metadata.calculation_performed, false);
  assert.equal(reference.formal_result_eligible, false); assert.equal(reference.predicted_synergy, null);
  assert.equal(reference.network_statistics.valid, false); assert.equal(reference.mechanism_type, 'not_computed');
  assert.match(G.modeLabel(reference), /非模型计算结果/); assert.equal(G.isSynthetic(reference), false);
});

test('every reference PPI is backed by a same-species UniProt INTERACTION record, not a pathway-derived edge', () => {
  const entries = reference.reference_provenance.uniprot_entries;
  const genes = new Map(entries.map(entry => [entry.accession, entry.gene]));
  const supported = new Set();
  for (const entry of entries) {
    assert.equal(entry.taxon_id, 9606); assert.match(entry.source_sha256, /^[a-f0-9]{64}$/);
    assert.equal(new URL(entry.source_url).hostname, 'rest.uniprot.org');
    for (const record of entry.selected_interactions) {
      const a = record.interactantOne.uniProtKBAccession, b = record.interactantTwo.uniProtKBAccession;
      assert.equal(a, entry.accession); assert.equal(record.organismDiffer, false);
      assert.ok(record.numberOfExperiments > 0); assert.ok(genes.has(a) && genes.has(b)); assert.notEqual(a, b);
      supported.add(pair('gene:' + genes.get(a), 'gene:' + genes.get(b)));
    }
  }
  const graph = reference.evidence_graphs.explanation_evidence_graph;
  const ppi = graph.edges.filter(edge => edge.type === 'PPI');
  assert.equal(ppi.length, 9); assert.deepEqual(new Set(ppi.map(edge => pair(edge.source, edge.target))), supported);
  for (const edge of ppi) {
    assert.equal(edge.directed, false); assert.equal(edge.confidence, null);
    for (const evidence of edge.evidence) {
      const source = entries.find(entry => entry.accession === evidence.entry);
      assert.ok(source); assert.equal(evidence.source_sha256, source.source_sha256);
      assert.ok(source.selected_interactions.some(record => JSON.stringify(record) === JSON.stringify(evidence.interaction_record)));
      assert.equal(pair(edge.source, edge.target), pair('gene:' + genes.get(evidence.interaction_record.interactantOne.uniProtKBAccession), 'gene:' + genes.get(evidence.interaction_record.interactantTwo.uniProtKBAccession)));
    }
    assert.ok(edge.evidence.length > 0);
  }
  const model = G.build(reference);
  assert.equal(model.nodeCount, 11); assert.equal(model.paths.length, 6); assert.ok(model.paths.every(row => row.valid));
});

test('reference memberships remain descriptive and illustrative positions cover existing nodes', () => {
  const graph = reference.evidence_graphs.explanation_evidence_graph;
  const members = reference.reference_provenance.reactome.selected_participants;
  for (const edge of graph.edges.filter(row => row.type === 'belongs_to_pathway')) assert.equal(members[edge.source.slice(5)], edge.reference_accession);
  const context = graph.edges.filter(row => row.type === 'belongs_to_disease_module');
  assert.equal(context.length, 1); assert.equal(context[0].source, 'gene:BRAF');
  assert.match(context[0].membership_basis, /single-gene literature background/);
  assert.deepEqual(new Set(Object.keys(reference.metadata.visualization_positions)), new Set(graph.nodes.map(node => node.id)));
  assert.ok(Object.values(reference.metadata.visualization_positions).every(point => point.length === 2 && point.every(Number.isFinite)));
  assert.ok(reference.metadata.licenses.some(source => source.license === 'CC BY 4.0'));
  assert.ok(reference.metadata.licenses.some(source => source.license === 'CC0'));
});
