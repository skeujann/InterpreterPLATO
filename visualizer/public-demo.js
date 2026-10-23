/* Synthetic UI demonstration, authored independently of local run files or databases. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./example.js'));
  else root.InterpreterPublicDemo = factory(root.InterpreterExample);
})(typeof globalThis !== 'undefined' ? globalThis : this, function (fixture) {
  'use strict';
  const payload = JSON.parse(JSON.stringify(fixture));
  payload.cell_line = 'DEMO';
  payload.run_mode = 'demo';
  payload.metadata.notice = '公开界面演示。全部节点和关系为合成数据，未运行 PLATO、解释器或统计检验。';
  payload.disease_context_source = 'synthetic_fixture_only';
  payload.disease_context_labels = ['synthetic_context'];
  payload.disease_context_size = 1;
  const graph = payload.evidence_graphs.explanation_evidence_graph;
  const added = Array.from({ length: 180 }, (_, i) => 'SYN_' + String(i + 1).padStart(3, '0'));
  graph.nodes.push(...added.map(name => ({ id: 'gene:' + name, type: 'gene', name, synthetic_fixture: true })));
  const pairs = new Set();
  function connect(left, right) {
    const key = [left, right].sort().join('|');
    if (left === right || pairs.has(key)) return;
    pairs.add(key);
    graph.edges.push({ source: 'gene:' + left, target: 'gene:' + right, type: 'PPI', directed: false,
      synthetic_fixture: true, source_dataset: 'synthetic_fixture_only_not_observed_ppi' });
  }
  // Fixed ring/chord pattern exercises overview and hover; no biological topology is used.
  added.forEach((name, i) => {
    connect(name, added[(i + 1) % added.length]);
    connect(name, added[(i + 13) % added.length]);
    if (i % 6 === 0) connect(name, i % 12 === 0 ? 'SYN_BRIDGE_1' : 'SYN_BRIDGE_2');
  });
  return {
    schema: 'interpreterplato-public-demo/v1',
    source: { case_id: 'PUBLIC DEMO', kind: 'synthetic', filename: 'public-demo.js', sha256: null },
    presentation: {
      nodes: ['drug:Demo-A', 'drug:Demo-B', 'gene:SYN_A', 'gene:SYN_B', 'gene:SYN_BRIDGE_1', 'gene:SYN_BRIDGE_2', 'gene:SYN_DISEASE', 'disease_module:SYN_CONTEXT'],
      pathIds: ['disease:0', 'disease:1'], endpoint: 'disease_module:SYN_CONTEXT',
      targets: ['gene:SYN_A', 'gene:SYN_B'], members: ['gene:SYN_DISEASE'],
      positions: { 'drug:Demo-A': [-430, -110], 'drug:Demo-B': [-430, 110],
        'gene:SYN_A': [-130, -110], 'gene:SYN_B': [-130, 110],
        'gene:SYN_BRIDGE_1': [160, -110], 'gene:SYN_BRIDGE_2': [160, 110],
        'gene:SYN_DISEASE': [450, 0], 'disease_module:SYN_CONTEXT': [810, 0] },
      compactPositions: { 'drug:Demo-A': [-120, -240], 'drug:Demo-B': [120, -240],
        'gene:SYN_A': [-120, -80], 'gene:SYN_B': [120, -80],
        'gene:SYN_BRIDGE_1': [-120, 70], 'gene:SYN_BRIDGE_2': [120, 70],
        'gene:SYN_DISEASE': [0, 220], 'disease_module:SYN_CONTEXT': [0, 390] },
      contextLabel: '合成疾病集合', contextSource: 'synthetic_fixture_only', contextSize: 1,
      contextDisplayLabel: '合成疾病集合\n非真实疾病注释',
    },
    payload,
  };
});
