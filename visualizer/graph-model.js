/* JSON-to-graph adapter. It never infers graph edges from paths or labels. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.InterpreterGraph = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const LIMITS = Object.freeze({ nodes: 5000, edges: 20000, results: 100, fileBytes: 20 * 1024 * 1024 });
  const COLORS = Object.freeze({ ppi: '#278b8b', target: '#c48638', disease: '#9468bd', pathway: '#ac94c6', report: '#a2aabc', hypothesis: '#c47587', other: '#8e99aa' });
  const NODE_TYPES = new Set(['drug', 'gene', 'pathway', 'disease_module']);
  const REPORT_TYPES = new Set(['supports_prediction', 'summarizes_combination_member', 'tests_sensitivity_of', 'changes_mechanism']);
  const TARGET_TYPES = new Set(['PLATO_predicted_target', 'known_target', 'known_and_predicted_target']);
  function object(value) { return value !== null && typeof value === 'object' && !Array.isArray(value); }
  function text(value, fallback) { return typeof value === 'string' && value.trim() ? value : fallback; }
  function relationship(type) {
    if (type === 'PPI') return 'ppi';
    if (TARGET_TYPES.has(type)) return 'target';
    if (type === 'belongs_to_disease_module') return 'disease';
    if (type === 'belongs_to_pathway') return 'pathway';
    if (REPORT_TYPES.has(type)) return 'report';
    if (type === 'hypothesized_priming_relation' || type === 'hypothesized_response_relation') return 'hypothesis';
    return 'other';
  }
  function records(input) {
    const rows = Array.isArray(input) ? input : (object(input) && Array.isArray(input.results) ? input.results : [input]);
    if (!rows.length || rows.length > LIMITS.results) throw new Error('解释记录必须为 1–100 条。');
    for (const row of rows) {
      if (!object(row) || !(row.error || row.evidence_graphs || row.result_schema_version || row.drug1_name)) {
        throw new Error('请选择解释结果 JSON；缓存、运行清单或任意网络 JSON 不属于此查看器的输入。');
      }
    }
    return rows;
  }
  function isSynthetic(payload) {
    return payload.generated_from_test_fixture === true ||
      (object(payload.metadata) && payload.metadata.generated_from_test_fixture === true);
  }
  function modeLabel(payload) {
    if (object(payload.metadata) && payload.metadata.conceptual_example === true) return '概念示意（非真实计算结果）· 合成测试数据';
    if (isSynthetic(payload)) return '合成测试数据 · 非真实机制证据';
    if (payload.run_mode === 'test' || payload.run_mode === 'demo') return '测试／演示模式 · 非正式机制结果';
    return payload.run_mode === 'research_hypothesis' ? '研究假说 · 来源未独立核验 · 非因果证据' : '运行模式未知 · 来源未独立核验 · 非因果证据';
  }
  function pairKey(left, right) { return JSON.stringify([left, right].sort()); }
  function build(payload) {
    if (!object(payload)) throw new Error('单条解释必须是 JSON 对象。');
    const graph = object(payload.evidence_graphs) && payload.evidence_graphs.explanation_evidence_graph;
    const warnings = [];
    const base = { payload, warnings, elements: [], paths: [], mode: modeLabel(payload), synthetic: isSynthetic(payload), roots: [], nodeCount: 0, ppiCount: 0 };
    if (!object(graph)) {
      warnings.push(payload.error ? '该记录未生成解释：' + String(payload.error) : 'JSON 不含 explanation_evidence_graph；不会从路径或靶点列表补画关系。');
      return base;
    }
    if (!Array.isArray(graph.nodes) || !Array.isArray(graph.edges)) throw new Error('证据图必须包含 nodes 和 edges 数组。');
    if (graph.nodes.length > LIMITS.nodes || graph.edges.length > LIMITS.edges) throw new Error('图超过本地查看上限：5,000 个节点／20,000 条边。');
    const nodes = new Map(), ppi = new Map(), elements = [], left = new Set(), right = new Set();
    for (const [side, targetSet] of [['drug1_targets', left], ['drug2_targets', right]]) {
      const targets = object(payload[side]) && payload[side].targets;
      if (Array.isArray(targets)) targets.forEach(target => { if (object(target) && typeof target.target_id === 'string') targetSet.add('gene:' + target.target_id); });
    }
    for (const [index, raw] of graph.nodes.entries()) {
      if (!object(raw) || typeof raw.id !== 'string' || !raw.id || raw.id.length > 2048) throw new Error('每个节点需要非空字符串 ID。');
      if (nodes.has(raw.id)) throw new Error('证据图含重复节点 ID：' + raw.id);
      const type = text(raw.type, 'other');
      const group = NODE_TYPES.has(type) ? type : 'report';
      const arm = left.has(raw.id) && right.has(raw.id) ? 'both' : left.has(raw.id) ? 'left' : right.has(raw.id) ? 'right' : 'bridge';
      const data = { id: 'n:' + index, originalId: raw.id, label: text(raw.name, raw.id).slice(0, 160), type, group, arm, raw };
      nodes.set(raw.id, data);
      elements.push({ group: 'nodes', data, classes: [group, 'arm-' + arm] });
      if (type === 'drug') base.roots.push(data.id);
    }
    let unknown = 0;
    for (const [index, raw] of graph.edges.entries()) {
      if (!object(raw) || !nodes.has(raw.source) || !nodes.has(raw.target)) throw new Error('证据图中的边引用了缺失节点；请检查上游 JSON。');
      const type = text(raw.type, 'unclassified');
      const category = relationship(type);
      const source = nodes.get(raw.source), target = nodes.get(raw.target);
      if (category === 'ppi' && (source.type !== 'gene' || target.type !== 'gene' || raw.directed === true || raw.source === raw.target)) {
        throw new Error('PPI 必须是两个不同基因节点之间的无向关系。');
      }
      const data = { id: 'e:' + index, source: source.id, target: target.id, type, category, raw };
      elements.push({ group: 'edges', data, classes: [category] });
      if (category === 'ppi') {
        const key = pairKey(raw.source, raw.target);
        if (!ppi.has(key)) ppi.set(key, []);
        ppi.get(key).push(data.id);
        base.ppiCount += 1;
      }
      if (category === 'other') unknown += 1;
    }
    if (unknown) warnings.push(unknown + ' 条未识别类型的边按未分类关系展示，不视为 PPI。');
    const local = object(payload.local_evidence) ? payload.local_evidence : {};
    const sources = [['cross', local.cross_drug_paths], ['disease', local.top_paths || payload.top_evidence_paths]];
    for (const [kind, rows] of sources) {
      if (!Array.isArray(rows)) continue;
      for (const [index, raw] of rows.entries()) {
        if (!object(raw) || !Array.isArray(raw.path) || !raw.path.length || raw.path.some(value => typeof value !== 'string' || !value)) {
          warnings.push('忽略一条缺少有效基因序列的路径记录。');
          continue;
        }
        const sequence = raw.path.map(gene => 'gene:' + gene);
        const missing = [], nodeIds = [], edgeIds = [];
        for (const originalId of sequence) {
          if (nodes.has(originalId)) nodeIds.push(nodes.get(originalId).id);
          else missing.push('缺少节点 ' + originalId);
        }
        for (let i = 1; i < sequence.length; i += 1) {
          const existing = ppi.get(pairKey(sequence[i - 1], sequence[i]));
          if (existing) edgeIds.push(...existing);
          else missing.push('缺少 PPI ' + raw.path[i - 1] + ' — ' + raw.path[i]);
        }
        base.paths.push({ id: kind + ':' + index, kind, raw, genes: raw.path.slice(), nodeIds, edgeIds: [...new Set(edgeIds)], missing,
          zeroHop: raw.path.length === 1, valid: missing.length === 0 });
      }
    }
    base.elements = elements;
    base.nodeCount = nodes.size;
    return base;
  }
  const styles = [
    { selector: 'node', style: { 'background-color': '#cad5df', 'border-color': '#ffffff', 'border-width': 2, label: 'data(label)', color: '#34465c', 'font-size': 12, 'font-family': 'Segoe UI, Microsoft YaHei, Arial, sans-serif', 'text-valign': 'bottom', 'text-margin-y': 9, 'text-wrap': 'wrap', 'text-max-width': 128, width: 32, height: 32 } },
    { selector: 'node.drug', style: { shape: 'round-rectangle', width: 116, height: 42, 'background-color': '#3b5b80', color: '#ffffff', 'text-valign': 'center', 'text-margin-y': 0, 'font-weight': 'bold', 'font-size': 12 } },
    { selector: 'node.gene.arm-left', style: { 'background-color': '#60a9a5' } },
    { selector: 'node.gene.arm-right', style: { 'background-color': '#dbad68' } },
    { selector: 'node.gene.arm-both', style: { 'background-color': '#a58ac4', 'border-color': '#725695', 'border-width': 3 } },
    { selector: 'node.disease_module', style: { shape: 'diamond', width: 55, height: 55, 'background-color': '#b496d2' } },
    { selector: 'node.pathway', style: { shape: 'hexagon', width: 43, height: 43, 'background-color': '#d5c5e8' } },
    { selector: 'node.report', style: { shape: 'round-rectangle', width: 30, height: 26, 'background-color': '#dbe0e8', 'border-style': 'dashed' } },
    { selector: 'edge', style: { width: 1.5, 'line-color': '#a2aabc', 'curve-style': 'bezier', 'source-arrow-shape': 'none', 'target-arrow-shape': 'none', 'mid-source-arrow-shape': 'none', 'mid-target-arrow-shape': 'none', opacity: 0.8 } },
    { selector: 'edge.ppi', style: { width: 2.3, 'line-color': COLORS.ppi, 'line-style': 'solid' } },
    { selector: 'edge.target', style: { 'line-color': COLORS.target, 'line-style': 'dashed' } },
    { selector: 'edge.disease', style: { 'line-color': COLORS.disease, 'line-style': 'dotted', width: 2 } },
    { selector: 'edge.pathway', style: { 'line-color': COLORS.pathway, 'line-style': 'dotted' } },
    { selector: 'edge.report', style: { 'line-color': COLORS.report, 'line-style': 'dashed', width: 1 } },
    { selector: 'edge.hypothesis', style: { 'line-color': COLORS.hypothesis, 'line-style': 'dashed' } },
    { selector: 'edge.other', style: { 'line-color': COLORS.other, 'line-style': 'dashed', opacity: 0.55 } },
    { selector: '.dimmed', style: { opacity: 0.13, 'text-opacity': 0.25 } },
    { selector: 'node.highlighted', style: { 'border-color': '#1a73e8', 'border-width': 4, 'z-index': 20 } },
    { selector: 'edge.highlighted', style: { 'line-color': '#1a73e8', width: 3.5, opacity: 1, 'z-index': 15 } },
    { selector: '.filtered', style: { display: 'none' } },
  ];
  return Object.freeze({ LIMITS, COLORS, records, build, styles, relationship, isSynthetic, modeLabel });
});
