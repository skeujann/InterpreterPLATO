/* Browser-only viewer. No fetch, upload endpoint, telemetry, or model execution. */
(function () {
  'use strict';
  const G = window.InterpreterGraph;
  const $ = id => document.getElementById(id);
  const RUNS = Array.isArray(window.InterpreterRunExamples) ? window.InterpreterRunExamples : window.InterpreterRunExample ? [window.InterpreterRunExample] : [];
  const QUALITY = { adequate: '充分', limited: '有限', insufficient: '不足' };
  const RULES = { network_close_distinct_processes: '网络接近，过程注释不同', network_proximal_targets: '候选靶点在 PPI 网络中接近', shared_pathway_co_modulation: '两侧候选靶点共同富集于部分通路', overlapping_single_exposure: '一侧疾病邻近，靶点或网络邻域存在重叠' };
  const state = { cy: null, rows: [], model: null, source: '', selectedPath: null, reference: false, runExample: false, referenceView: 'focus', compact: false };
  const COMPACT_POSITIONS = { 'drug:Dabrafenib': [-100, -220], 'drug:Trametinib': [100, -220], 'gene:MAP2K1': [40, -80], 'gene:MAP2K2': [150, 10], 'gene:BRAF': [-65, 75], 'disease_module:literature-context': [50, 240] };
  const GROUP_LABELS = { gene: '基因／蛋白节点', drug: '药物', pathway: '通路集合', disease_module: '疾病基因集合', prediction_evidence: '预测证据记录', mechanism: '报告摘要结构', counterfactual: '敏感性分析记录' };
  const EDGE_LABELS = { ppi: '文件记录的无向 PPI', target: '候选靶点关联', disease: '疾病集合归属', pathway: '通路集合归属', report: '报告结构／证据支持', hypothesis: '规则生成假说', other: '未分类关系' };
  function element(tag, className, value) {
    const result = document.createElement(tag);
    if (className) result.className = className;
    if (value !== undefined) result.textContent = String(value);
    return result;
  }
  function error(message) { $('error-banner').textContent = message; $('error-banner').hidden = !message; }
  function download(blob, name) {
    const url = URL.createObjectURL(blob), link = element('a');
    link.href = url; link.download = name; document.body.appendChild(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function pretty(value, limit) {
    const result = JSON.stringify(value, null, 2);
    return result.length > limit ? result.slice(0, limit) + '\n…（界面截断；完整属性保留在输入文件）' : result;
  }
  function activateTab(name, names) {
    for (const target of names) {
      $('tab-' + target).setAttribute('aria-selected', String(name === target));
      $('tab-' + target).tabIndex = name === target ? 0 : -1;
      $(target + '-panel').hidden = name !== target;
    }
  }
  function tab(name) { activateTab(name, ['details', 'diagnostics']); }
  function exploreTab(name) { activateTab(name, ['paths', 'filters']); }
  function bindTabs(names, activate) {
    names.forEach((name, index) => {
      const button = $('tab-' + name);
      button.addEventListener('click', () => activate(name));
      button.addEventListener('keydown', event => {
        const target = event.key === 'ArrowRight' ? (index + 1) % names.length
          : event.key === 'ArrowLeft' ? (index + names.length - 1) % names.length
          : event.key === 'Home' ? 0 : event.key === 'End' ? names.length - 1 : -1;
        if (target < 0) return;
        event.preventDefault(); activate(names[target]); $('tab-' + names[target]).focus();
      });
    });
  }
  function inspector(open) {
    $('inspector-panel').hidden = !open;
    $('explorer-body').classList.toggle('inspector-open', open);
    $('toggle-inspector').setAttribute('aria-expanded', String(open));
    $('toggle-inspector').classList.toggle('active', open);
    if (state.cy) state.cy.resize();
  }
  function sidebar(open) {
    $('explore-panel').hidden = !open;
    $('explorer-body').classList.toggle('sidebar-collapsed', !open);
    $('toggle-sidebar').setAttribute('aria-expanded', String(open));
    $('toggle-sidebar').classList.toggle('active', open);
    if (state.cy) state.cy.resize();
  }
  function zoom(factor) {
    if (!state.cy) return;
    state.cy.zoom({ level: Math.max(state.cy.minZoom(), Math.min(state.cy.maxZoom(), state.cy.zoom() * factor)),
      renderedPosition: { x: $('cy').clientWidth / 2, y: $('cy').clientHeight / 2 } });
  }
  function referenceMode() { return state.model && state.model.payload.run_mode === 'reference_example'; }
  function shown(item) { return !item.hasClass('filtered') && !item.hasClass('focus-hidden'); }
  function guided() { return state.reference || state.runExample; }
  function focusView() { return guided() && state.referenceView === 'focus'; }
  function guide() {
    return state.runExample ? window.InterpreterRunExample.presentation : {
      nodes: window.InterpreterReference.metadata.focus_view.nodes,
      positions: window.InterpreterReference.metadata.visualization_positions, compactPositions: COMPACT_POSITIONS,
    };
  }
  function compactView() { return focusView() && $('cy').clientWidth > 0 && $('cy').clientWidth < 560; }
  function caption() {
    return focusView() ? '终点为疾病模块 · 紫色点线为' + (state.runExample ? '集合归属' : '背景归属') + ' · 连线无方向' : '仅展示文件中的关系 · 连线无方向';
  }
  function renderRunSummary() {
    $('run-summary').hidden = !state.runExample; $('run-guide').hidden = !state.runExample;
    $('run-suite').hidden = !state.runExample || RUNS.length < 2;
    document.querySelectorAll('[data-run-index]').forEach(button => button.setAttribute('aria-pressed', String(state.runExample && RUNS[Number(button.dataset.runIndex)].payload === state.model.payload)));
    if (!state.runExample) return;
    const p = state.model.payload, stats = p.network_statistics || {}, source = window.InterpreterRunExample.source;
    const evaluated = stats.valid === true;
    $('run-verdict').textContent = !evaluated ? '网络统计不可评估' : stats.statistically_supported === true
      ? typeof stats.near_A === 'boolean' && typeof stats.near_B === 'boolean' && stats.near_A !== stats.near_B ? '仅一侧疾病邻近获支持' : '网络统计记录：获得支持'
      : stats.statistically_supported === false ? '疾病邻近未获显著支持' : '统计支持状态未记录';
    $('run-rule').textContent = '探索性规则分类：' + (RULES[p.mechanism_type] || p.mechanism_type) + '。' + (stats.n_randomizations || '未记录') + ' 次随机化 · 统计质量：' + (QUALITY[stats.quality_status] || stats.quality_status || '未记录') + '。';
    const metrics = $('run-metrics'); metrics.replaceChildren();
    for (const [label, value, note] of [
      [p.drug1_name + ' / 校正 p', evaluated ? stats.p_A_adjusted : null, evaluated ? '疾病邻近检验' : '统计质量不足，未提供有效 p 值'],
      [p.drug2_name + ' / 校正 p', evaluated ? stats.p_B_adjusted : null, evaluated ? '疾病邻近检验' : '统计质量不足，未提供有效 p 值'],
      ['上游输入分数', p.predicted_synergy, '由调用方提供，非本次解释器计算'],
    ]) {
      const item = element('div'); item.append(element('dt', '', label), element('dd', '', typeof value === 'number' && Number.isFinite(value) ? value.toFixed(4) : !evaluated && label.includes('校正 p') ? '不可评估' : '未记录'), element('small', '', note)); metrics.appendChild(item);
    }
    $('run-guide-id').textContent = source.case_id;
    $('run-guide-note').textContent = '来自已有运行 JSON · ' + guide().pathIds.length + ' 条选定疾病路径';
  }
  function showRunDetails() {
    if (!state.runExample) return;
    const p = state.model.payload, s = p.network_statistics || {};
    $('detail-title').textContent = window.InterpreterRunExample.source.case_id + ' · 实际运行输出';
    $('detail-description').textContent = '开发测试结果；规则分类与统计支持分别呈现。疾病集合为输入中的 ' + guide().contextLabel + '，不表示细胞系特异疾病基因集。';
    $('detail-json').textContent = pretty({ source_file: window.InterpreterRunExample.source,
      run_mode: p.run_mode, generated_from_test_fixture: p.generated_from_test_fixture, formal_result_eligible: p.formal_result_eligible,
      mechanism_type: p.mechanism_type, predicted_synergy_upstream_input: p.predicted_synergy,
      disease_context: { source: p.disease_context_source, labels: p.disease_context_labels, size: p.disease_context_size },
      statistics: { method: s.method_version, statistically_supported: s.statistically_supported, valid: s.valid, quality_status: s.quality_status,
        randomizations: s.n_randomizations, raw_p_A: s.p_A, raw_p_B: s.p_B, adjusted_p_A: s.p_A_adjusted, adjusted_p_B: s.p_B_adjusted, correction: s.multiple_testing_method },
      model_fidelity: p.model_fidelity,
    }, 14000);
    $('detail-json').hidden = false; tab('details'); inspector(true);
  }
  function viewCopy() {
    $('reference-view-switch').hidden = !guided();
    $('reference-view-switch').setAttribute('aria-label', state.runExample ? '运行实例显示范围' : '文献实例显示范围');
    $('view-focus').setAttribute('aria-pressed', String(focusView()));
    $('view-full').setAttribute('aria-pressed', String(guided() && !focusView()));
    $('arm-legend').hidden = guided(); $('reference-legend').hidden = !state.reference; $('run-legend').hidden = !state.runExample;
    $('focus-stage-labels').hidden = !focusView() || compactView() || $('layout-select').value !== 'preset';
    $('path-section-title').textContent = focusView() ? '读到疾病模块' : '已记录路径';
    $('view-title').replaceChildren();
    $('view-title').appendChild(document.createTextNode(state.reference ? '从双药靶点，' : '探索连接，'));
    $('view-title').appendChild(element('br'));
    $('view-title').appendChild(document.createTextNode(state.reference ? '读到疾病背景。' : '读懂已有证据。'));
    $('view-description').textContent = state.reference ? '沿已有互作，经 BRAF 读到疾病模块的文献背景。' : '打开解释 JSON，在本地探索药物、靶点与证据网络。';
    $('hero-type').textContent = state.reference ? '阅读终点 / 疾病模块' : 'LOCAL EVIDENCE';
    $('hero-endpoint').replaceChildren(document.createTextNode(state.reference ? '黑色素瘤' : 'JSON'), element('span', 'endpoint-period', '.'));
    $('hero-endpoint').classList.toggle('disease-title', state.reference);
    $('hero-endpoint').classList.toggle('run-title', state.runExample);
    $('hero-note').textContent = state.reference ? 'BRAF V600 · 文献背景，非计算疾病基因集。' : '文件留在本地，关系来自输入记录。';
    $('canvas-caption-text').textContent = caption();
    const stages = state.runExample ? ['药物', 'PLATO 候选靶点', '疾病集合成员', '疾病模块'] : ['药物', '靶点与互作', 'BRAF', '疾病模块'];
    [...$('focus-stage-labels').children].forEach((span, index) => span.replaceChildren(element('b', '', String(index + 1).padStart(2, '0')), document.createTextNode(' ' + stages[index])));
    if (state.runExample) {
      $('view-title').replaceChildren(document.createTextNode('从 PLATO 候选，'), element('br'), document.createTextNode('读到真实运行结果。'));
      $('view-description').textContent = '读取已有解释器输出，分开展示网络路径、规则分类与统计结果。';
      $('hero-type').textContent = '分析细胞系 / 开发测试';
      $('hero-endpoint').replaceChildren(document.createTextNode(state.model.payload.cell_line || 'LOCAL'), element('span', 'endpoint-period', '.'));
      $('hero-endpoint').classList.remove('disease-title');
      $('hero-note').textContent = guide().contextLabel + ' · ' + guide().contextSize + ' 个基因 · ' + window.InterpreterRunExample.source.case_id;
    }
    renderRunSummary();
  }
  function setReferenceView(kind) {
    if (!guided() || !state.cy) return;
    clearHighlight(); state.referenceView = kind === 'full' ? 'full' : 'focus';
    resetFilters(focusView());
    if (focusView()) {
      document.querySelector('[data-node-group="disease_module"]').checked = true;
      document.querySelector('[data-edge-group="membership"]').checked = true;
    }
    if (state.runExample && !focusView()) document.querySelectorAll('[data-node-group], [data-edge-group]').forEach(input => { input.checked = true; });
    const ids = new Set(guide().nodes);
    state.cy.batch(() => {
      state.cy.nodes().forEach(node => node.toggleClass('focus-hidden', focusView() && !ids.has(node.data('originalId'))));
      state.cy.edges().forEach(edge => edge.toggleClass('focus-hidden', focusView() &&
        (edge.source().hasClass('focus-hidden') || edge.target().hasClass('focus-hidden') || !['ppi', 'target', 'disease'].includes(edge.data('category')))));
    });
    $('node-search').value = ''; $('search-status').textContent = '输入名称可高亮匹配节点';
    $('layout-select').value = state.runExample && !focusView() ? 'cose' : 'preset'; viewCopy(); applyFilters(); renderPaths(); inspector(false); layout();
  }
  function edgeLabel(category) {
    return referenceMode() && category === 'disease' ? '文献背景归属'
      : referenceMode() && category === 'target' ? '文献靶点关联' : EDGE_LABELS[category];
  }
  function inspect(item) {
    const data = item.data();
    $('detail-title').textContent = item.isNode() ? data.label : edgeLabel(data.category);
    $('detail-description').textContent = item.isNode()
      ? (referenceMode() && data.type === 'disease_module' ? '文献背景注释（范围见原始属性）' : GROUP_LABELS[data.type] || data.type) + ' · ' + data.originalId
      : data.type + (data.category === 'ppi' ? ' · 无方向；来源未在此独立核验' : ' · 关联或结构关系，不作因果箭头解读');
    $('detail-json').textContent = pretty(data.raw, 6500);
    $('detail-json').hidden = false;
  }
  function clearHighlight() {
    state.selectedPath = null;
    if (state.cy) state.cy.elements().removeClass('dimmed highlighted');
    document.querySelectorAll('.path-button').forEach(button => { button.classList.remove('selected'); button.setAttribute('aria-pressed', 'false'); });
    document.querySelectorAll('[data-tour]').forEach(button => button.setAttribute('aria-pressed', 'false'));
  }
  function resetFilters(reference) {
    document.querySelectorAll('[data-node-group]').forEach(input => {
      input.checked = ['drug', 'gene'].includes(input.dataset.nodeGroup) || (!reference && input.dataset.nodeGroup !== 'report');
    });
    document.querySelectorAll('[data-edge-group]').forEach(input => {
      input.checked = ['ppi', 'target'].includes(input.dataset.edgeGroup) || (!reference && input.dataset.edgeGroup === 'membership');
    });
  }
  function builtin() {
    const selected = $('example-select').value;
    const index = selected === 'run' ? 0 : selected.startsWith('run:') ? Number(selected.slice(4)) : -1;
    if (Number.isInteger(index) && RUNS[index]) { window.InterpreterRunExample = RUNS[index]; return RUNS[index].payload; }
    return selected === 'synthetic' ? window.InterpreterExample : window.InterpreterReference;
  }
  function loadBuiltin() {
    const payload = builtin(), reference = payload === window.InterpreterReference;
    $('case-guide').open = false;
    resetFilters(reference); sidebar(true); exploreTab('paths');
    const actualRun = window.InterpreterRunExample && payload === window.InterpreterRunExample.payload;
    load(payload, actualRun ? window.InterpreterRunExample.source.case_id + ' · SHA-256 ' + window.InterpreterRunExample.source.sha256.slice(0, 12) : reference ? '文献实例 · 公共注释快照 2026-10-21' : '内置合成概念示例 · 无真实计算');
    if (reference) tour('overview');
  }
  function referenceSources() {
    const container = $('reference-sources'); container.replaceChildren();
    container.appendChild(element('h3', '', '来源与整理范围'));
    container.appendChild(element('p', 'reference-limit', '人工整理的文献与数据库注释。PPI 不一定是直接结合，也未在该突变或黑色素瘤背景逐边验证。BRAF 疾病背景是单基因文献注释，不是疾病基因集或富集结果。'));
    for (const source of window.InterpreterReference.metadata.references) {
      const card = element('div', 'source-card'), link = element('a', '', source.title + ' ↗');
      link.href = source.url; link.target = '_blank'; link.rel = 'noopener noreferrer';
      card.appendChild(link); card.appendChild(element('p', '', source.detail)); container.appendChild(card);
    }
    for (const source of window.InterpreterReference.metadata.licenses) {
      const item = element('p', 'reference-license', source.resource + ' · ');
      const link = element('a', '', source.license); link.href = source.url; link.target = '_blank'; link.rel = 'noopener noreferrer';
      item.appendChild(link); container.appendChild(item);
    }
  }
  function showSources() {
    if (!state.reference) return;
    $('detail-title').textContent = '文献实例的来源';
    $('detail-description').textContent = '真实药物与公开注释；没有运行 PLATO、解释器或网络统计。';
    $('detail-json').hidden = true; tab('details'); inspector(true);
  }
  function tour(kind) {
    if (!state.reference || !state.cy) return;
    setReferenceView(kind === 'context' ? 'full' : 'focus');
    if (kind === 'overview') {
      $('tour-status').textContent = '先看双药与 MEK1 / MEK2，再经 BRAF 读到最右侧疾病模块。紫色点线仅为文献背景归属。';
      inspector(false);
    } else if (kind === 'ppi') {
      const edges = state.cy.edges().filter(edge => shown(edge) && edge.data('category') === 'ppi');
      state.cy.elements().addClass('dimmed'); edges.union(edges.connectedNodes()).removeClass('dimmed').addClass('highlighted');
      $('detail-title').textContent = 'MEK1 / MEK2 — BRAF';
      $('detail-description').textContent = '2 条公开记录的无向 PPI 连接 MEK1 / MEK2 与 BRAF；BRAF 另由背景归属点线关联疾病模块。不表示激活、抑制或信号方向。';
      $('detail-json').hidden = true; tab('details'); inspector(true);
      $('tour-status').textContent = '2 条已有 PPI，连接 MEK1 / MEK2 与 BRAF。是网络连接，不是因果或协同结论。';
    } else {
      const edges = state.cy.edges().filter(edge => kind === 'targets' ? edge.data('category') === 'target' : ['disease', 'pathway'].includes(edge.data('category')));
      state.cy.elements().addClass('dimmed'); edges.union(edges.connectedNodes()).removeClass('dimmed').addClass('highlighted');
      $('detail-title').textContent = kind === 'targets' ? '两个药物靶点层级' : '通路与文献背景';
      $('detail-description').textContent = kind === 'targets'
        ? 'Dabrafenib 关联 BRAF；Trametinib 关联 MEK1 / MEK2。这里只展示选定的药理靶点，不是 PLATO 预测。'
        : '紫色点线表示集合归属。BRAF V600 黑色素瘤为文献背景单基因注释，未进行疾病或通路富集。';
      $('detail-json').hidden = true; tab('details'); inspector(true);
      $('tour-status').textContent = kind === 'targets' ? '3 条有出处的药物—靶点关联；并不表示完整靶点谱或细胞系实测。' : '通路与疾病背景只作注释展示，不能与 PPI 混为作用箭头。';
    }
    document.querySelectorAll('[data-tour]').forEach(button => button.setAttribute('aria-pressed', String(button.dataset.tour === kind)));
  }
  function renderCounts() {
    if (!state.cy) return;
    const visible = state.cy.nodes().filter(shown);
    const ppi = state.cy.edges().filter(edge => edge.data('category') === 'ppi' && shown(edge));
    $('graph-counts').textContent = visible.length + ' / ' + state.cy.nodes().length + ' 节点 · ' + ppi.length + ' 条 PPI';
    $('graph-empty').hidden = visible.length > 0;
    $('graph-empty').querySelector('strong').textContent = state.cy.nodes().length ? '当前筛选未显示任何节点' : '该记录没有可展示的证据图';
    $('export-png').disabled = !visible.length;
    for (const id of ['zoom-in', 'zoom-out', 'fit-graph']) $(id).disabled = !visible.length;
  }
  function applyFilters() {
    if (!state.cy) return;
    const nodeGroups = new Set([...document.querySelectorAll('[data-node-group]:checked')].map(input => input.dataset.nodeGroup));
    const edgeGroups = new Set([...document.querySelectorAll('[data-edge-group]:checked')].map(input => input.dataset.edgeGroup));
    state.cy.batch(() => {
      state.cy.nodes().forEach(node => node.toggleClass('filtered', !nodeGroups.has(node.data('group'))));
      state.cy.edges().forEach(edge => {
        const category = edge.data('category');
        const group = ['disease', 'pathway'].includes(category) ? 'membership' : ['report', 'hypothesis', 'other'].includes(category) ? 'report' : category;
        edge.toggleClass('filtered', !edgeGroups.has(group) || edge.source().hasClass('filtered') || edge.target().hasClass('filtered'));
      });
    });
    renderCounts();
  }
  function layout() {
    if (!state.cy || !state.model.nodeCount) return;
    const name = $('layout-select').value;
    state.compact = compactView();
    state.cy.nodes().toggleClass('compact-reference', state.compact);
    $('focus-stage-labels').hidden = !focusView() || state.compact || name !== 'preset';
    $('canvas-caption-text').textContent = caption();
    const options = { name, animate: false, fit: true, padding: 40, nodeDimensionsIncludeLabels: true };
    if (name === 'preset') {
      if (!guided()) { $('layout-select').value = 'cose'; layout(); return; }
      const positions = state.compact ? guide().compactPositions : guide().positions;
      state.cy.nodes().forEach(node => { const point = positions[node.data('originalId')]; if (point) node.position({ x: point[0], y: point[1] }); });
      options.padding = state.compact ? 28 : focusView() ? 45 : 48;
    }
    const visibleNodes = state.cy.nodes().filter(shown).length;
    if (name === 'cose' && visibleNodes > 400) {
      options.name = 'grid'; $('layout-select').value = 'grid';
      $('search-status').textContent = '节点较多，使用轻量网格布局。可筛选后重新布局。';
    }
    if (options.name === 'cose') Object.assign(options, { randomize: false, nodeRepulsion: 16000, idealEdgeLength: 110, numIter: 800 });
    if (name === 'breadthfirst') Object.assign(options, { directed: false, roots: state.cy.nodes().filter(node => node.data('type') === 'drug'), spacingFactor: 1.4 });
    if (name === 'concentric') Object.assign(options, { minNodeSpacing: 35, concentric: node => node.data('type') === 'drug' ? 10 : node.degree(), levelWidth: () => 2 });
    state.cy.elements().filter(shown).layout(options).run();
  }
  function selectPath(path, button) {
    clearHighlight(); $('node-search').value = '';
    if (!path.valid) return;
    state.selectedPath = path.id;
    const ids = new Set([...path.nodeIds, ...path.edgeIds]);
    // The curated reading route may include its actual drug-target edge. Never infer one.
    if (guided() && path.raw.drug) {
      const drug = state.cy.nodes().filter(node => node.data('type') === 'drug' && node.data('label') === path.raw.drug);
      const first = path.nodeIds[0];
      const links = state.cy.edges().filter(edge => edge.data('category') === 'target' &&
        ((drug.contains(edge.source()) && edge.target().id() === first) || (drug.contains(edge.target()) && edge.source().id() === first)));
      links.union(links.connectedNodes()).forEach(item => ids.add(item.id()));
    }
    if (guided() && path.kind === 'disease') {
      const last = path.nodeIds[path.nodeIds.length - 1];
      const memberships = state.cy.edges().filter(edge => edge.data('category') === 'disease' &&
        ((edge.source().id() === last && edge.target().data('type') === 'disease_module') ||
         (edge.target().id() === last && edge.source().data('type') === 'disease_module')));
      memberships.union(memberships.connectedNodes()).forEach(item => ids.add(item.id()));
    }
    state.cy.elements().addClass('dimmed');
    const selected = state.cy.elements().filter(item => ids.has(item.id()));
    selected.removeClass('dimmed').addClass('highlighted');
    button.classList.add('selected'); button.setAttribute('aria-pressed', 'true');
    const hidden = selected.filter(item => !shown(item)).length;
    $('search-status').textContent = path.zeroHop ? '零跳记录：没有 PPI 边，也不会补画自环。' : '高亮 ' + path.edgeIds.length + ' 条文件中已有的 PPI 边。';
    if (hidden) $('search-status').textContent += ' 当前筛选隐藏了部分路径元素。';
    $('detail-title').textContent = focusView() ? (path.raw.drug || '靶点') + ' / 疾病模块' : path.kind === 'cross' ? '双药靶点连接' : referenceMode() ? '靶点—文献背景' : '靶点—疾病路径';
    $('detail-description').textContent = focusView() ? state.runExample ? '高亮原始运行中的候选靶点关联、无向 PPI 与疾病集合归属；没有增添边，也不将连通解释为统计显著性。' : '已有基因路径与疾病背景归属分开呈现。紫色点线不是 PPI；该背景不是计算疾病基因集。' : path.zeroHop ? '零跳仅表示共享节点或集合成员关系。' : '路径仅用于高亮已记录的边，不增添网络关系。';
    $('detail-json').textContent = pretty(path.raw, 6500); $('detail-json').hidden = false; tab('details'); inspector(true);
  }
  function renderPaths() {
    const list = $('path-list'); list.replaceChildren();
    const displayed = state.model.paths.filter(path => !focusView() || (state.runExample ? guide().pathIds.includes(path.id) : path.kind === 'disease' && path.genes[path.genes.length - 1] === 'BRAF'));
    $('path-count').textContent = String(displayed.length);
    if (!displayed.length) { list.appendChild(element('p', 'small muted', '该记录未包含可展示的路径。')); return; }
    for (const [kind, title] of [['cross', '双药靶点连接'], ['disease', referenceMode() ? '靶点—文献背景' : '靶点—疾病路径']]) {
      const paths = displayed.filter(path => path.kind === kind);
      if (!paths.length) continue;
      list.appendChild(element('div', 'path-group-label', focusView() ? 'READING ROUTES / ' + String(displayed.length).padStart(2, '0') : title + ' · ' + paths.length));
      for (const [index, path] of paths.entries()) {
        const button = element('button', 'path-button'); button.type = 'button'; button.disabled = !path.valid;
        button.dataset.pathId = path.id;
        button.setAttribute('aria-pressed', 'false');
        const heading = element('span', 'path-heading');
        heading.appendChild(element('span', 'path-number', String(index + 1).padStart(2, '0')));
        heading.appendChild(element('strong', '', path.raw.drug || '路径 ' + (index + 1)));
        button.appendChild(heading);
        button.appendChild(element('span', 'path-route', path.genes.map(gene => state.reference && gene === 'MAP2K1' ? 'MEK1' : state.reference && gene === 'MAP2K2' ? 'MEK2' : gene).join(' — ')));
        if (focusView()) button.appendChild(element('span', 'path-destination', state.runExample ? '⋯ ' + guide().contextLabel + ' · 疾病集合' : '⋯ 黑色素瘤 · 文献背景'));
        if (focusView()) button.appendChild(element('span', 'path-zero', path.zeroHop ? state.runExample ? '候选靶点即集合成员 · 0 条 PPI' : '直接靶点关联 · 0 条 PPI' : '靶点关联 + ' + path.edgeIds.length + ' 条无向 PPI'));
        else if (path.zeroHop) button.appendChild(element('span', 'path-zero', kind === 'cross' ? '共享节点零跳 · 0 条 PPI' : referenceMode() ? '文献背景成员零跳 · 0 条 PPI' : '疾病集合成员零跳 · 0 条 PPI'));
        else button.appendChild(element('span', 'path-zero', (path.genes.length - 1) + ' 跳 · 仅高亮文件中已有的边'));
        if (!path.valid) button.appendChild(element('span', 'path-issue', '无法高亮：' + path.missing.join('；')));
        button.addEventListener('click', () => selectPath(path, button)); list.appendChild(button);
      }
    }
  }
  function diagnostics() {
    const container = $('diagnostic-summary'); container.replaceChildren();
    for (const warning of state.model.warnings) container.appendChild(element('p', 'diagnostic-item', warning));
    const invalid = state.model.paths.filter(path => !path.valid).length;
    if (invalid) container.appendChild(element('p', 'diagnostic-item', invalid + ' 条路径与证据图不一致，已禁用高亮；没有自动补边。'));
    const payload = state.model.payload;
    if (Array.isArray(payload.warning_flags) && payload.warning_flags.length) container.appendChild(element('p', 'diagnostic-item', '上游记录了 ' + payload.warning_flags.length + ' 条提示，完整标记见下方 JSON。'));
    if (payload.error) container.appendChild(element('p', 'diagnostic-item', '上游错误：' + String(payload.error)));
    $('diagnostic-count').textContent = String(container.children.length);
    $('diagnostic-count').hidden = !container.children.length;
    $('diagnostic-json').textContent = pretty({
      network_connectivity: payload.network_connectivity || (payload.local_evidence && payload.local_evidence.connectivity) || { status: 'not_recorded' },
      warning_flags: payload.warning_flags || [],
      network_statistics_status: payload.network_statistics ? { valid: payload.network_statistics.valid, quality_status: payload.network_statistics.quality_status } : 'not_recorded',
      viewer_warnings: state.model.warnings,
    }, 16000);
  }
  function display(index, suppliedModel) {
    const payload = state.rows[index], model = suppliedModel || G.build(payload);
    state.model = model; state.selectedPath = null;
    state.reference = payload === window.InterpreterReference;
    state.runExample = !!window.InterpreterRunExample && payload === window.InterpreterRunExample.payload;
    state.referenceView = 'focus';
    $('case-guide').hidden = !state.reference; $('reference-sources').hidden = !state.reference;
    $('explore-panel').classList.toggle('has-guide', state.reference);
    $('reference-layout').hidden = !guided(); $('reference-layout').disabled = !guided();
    if (guided()) $('layout-select').value = 'preset';
    else if ($('layout-select').value === 'preset') $('layout-select').value = 'cose';
    if (state.cy) state.cy.destroy();
    state.cy = null; inspector(false); tab('details');
    $('node-search').value = ''; $('search-status').textContent = '输入名称可高亮匹配节点';
    $('case-title').textContent = (payload.drug1_name || payload.drug1 || '未标注药物 A') + ' + ' + (payload.drug2_name || payload.drug2 || '未标注药物 B');
    $('case-meta').textContent = state.reference ? '' : (payload.cell_line ? payload.cell_line + ' · ' : '') + '记录 ' + (index + 1) + ' / ' + state.rows.length;
    $('case-meta').hidden = state.reference;
    $('source-name').textContent = state.source;
    const concept = payload.metadata && payload.metadata.conceptual_example === true;
    $('mode-banner').className = 'mode-banner ' + (model.synthetic || payload.run_mode === 'test' || payload.run_mode === 'demo' ? 'synthetic' : 'research');
    $('mode-title').textContent = concept ? '概念示意（非真实计算结果）' : model.mode;
    $('mode-description').textContent = concept ? '全部节点与关系为合成测试数据，不用于科研证据。' : '只展示输入文件记录的关系，不独立核验原始 PPI，不推导因果或机制真实性。';
    if (payload.run_mode === 'reference_example') $('mode-description').textContent = '人工整理的药物与公共注释；互作未在该临床背景逐边验证，不运行模型或统计，不证明协同。';
    if (state.runExample) {
      $('mode-title').textContent = '实际运行实例 · 开发测试';
      $('mode-description').textContent = '读取已有解释器 JSON · PLATO 候选靶点 · 保留原始统计与诊断，非正式验证结果。';
    }
    state.cy = window.cytoscape({ container: $('cy'), elements: model.elements, style: G.styles, layout: { name: 'grid', padding: 40 }, minZoom: 0.12, maxZoom: 3.5 });
    if (state.reference) {
      state.cy.nodes().forEach(node => {
        const id = node.data('originalId');
        if (node.data('type') === 'drug') node.addClass('reference-drug');
        if (id === 'gene:BRAF') node.addClass('reference-anchor');
        if (id === 'disease_module:literature-context') {
          node.addClass('reference-context'); node.data('displayLabel', '疾病模块 / 黑色素瘤\nBRAF V600 · 文献背景');
        }
        if (['gene:MAP2K1', 'gene:MAP2K2'].includes(id)) {
          node.addClass('reference-target'); node.data('displayLabel', id === 'gene:MAP2K1' ? 'MEK1\nMAP2K1' : 'MEK2\nMAP2K2');
        }
      });
    }
    if (state.runExample) {
      const view = guide();
      state.cy.nodes().forEach(node => {
        const id = node.data('originalId');
        if (node.data('type') === 'drug') node.addClass('run-drug');
        if (view.targets.includes(id)) { node.addClass('run-target'); node.data('displayLabel', id.startsWith('gene:') ? id.slice(5) : node.data('label')); }
        if (view.members.includes(id)) { node.addClass('run-member'); node.data('displayLabel', id.startsWith('gene:') ? id.slice(5) : node.data('label')); }
        if (id === view.endpoint) { node.addClass('run-context'); node.data('displayLabel', view.contextDisplayLabel); }
      });
    }
    state.cy.on('tap', 'node, edge', event => { inspect(event.target); tab('details'); inspector(true); });
    state.cy.on('mouseover', 'node, edge', event => {
      const item = event.target; inspect(item);
      const position = event.renderedPosition || (item.isNode() ? item.renderedPosition() : item.renderedMidpoint());
      $('graph-tooltip').textContent = item.isNode() ? item.data('label') + '\n' + (referenceMode() && item.data('type') === 'disease_module' ? '文献背景注释' : GROUP_LABELS[item.data('type')] || item.data('type')) : edgeLabel(item.data('category')) + '\n' + item.data('type');
      $('graph-tooltip').style.left = Math.max(8, Math.min(position.x + 16, $('cy').clientWidth - 240)) + 'px';
      $('graph-tooltip').style.top = Math.max(8, Math.min(position.y + 18, $('cy').clientHeight - 90)) + 'px';
      $('graph-tooltip').hidden = false;
    });
    state.cy.on('mouseout', 'node, edge', () => { $('graph-tooltip').hidden = true; });
    $('detail-title').textContent = '点选一个节点'; $('detail-description').textContent = '点击节点、连线或左侧路径，查看文件中记录的属性。';
    $('detail-json').hidden = true; $('graph-tooltip').hidden = true;
    viewCopy(); applyFilters(); renderPaths(); diagnostics(); layout();
    if (guided()) setReferenceView('focus');
    if (!model.nodeCount) { tab('diagnostics'); inspector(true); }
  }
  function load(input, name) {
    const rows = G.records(input); G.build(rows[0]);
    state.rows = rows; state.source = name;
    $('result-select').replaceChildren();
    rows.forEach((row, index) => $('result-select').appendChild(new Option((index + 1) + '. ' + (row.drug1_name || row.drug1 || '药物 A') + ' + ' + (row.drug2_name || row.drug2 || '药物 B') + (row.error ? '（未完成）' : ''), String(index))));
    $('result-select').hidden = rows.length < 2; $('result-label').hidden = rows.length < 2;
    error(''); display(0);
  }
  async function openFile(file) {
    if (!file) return;
    try {
      if (file.size > G.LIMITS.fileBytes) throw new Error('文件超过 20 MiB。本查看器用于轻量解释结果，不读取原始数据集。');
      load(JSON.parse((await file.text()).replace(/^\uFEFF/, '')), file.name);
    } catch (exception) { error('文件未载入：' + exception.message + ' 当前视图仍为之前的记录。'); }
    $('file-input').value = '';
  }
  async function exportPng() {
    if (!state.cy || !state.model.nodeCount) return;
    const exportedModel = state.model, exportedCaseId = state.runExample ? window.InterpreterRunExample.source.case_id : null;
    const dataUrl = state.cy.png({ full: true, scale: 2, bg: '#fbfcf6', maxWidth: 2400, maxHeight: 1700 });
    const image = new Image();
    image.onload = () => {
      const canvas = document.createElement('canvas'); canvas.width = Math.max(image.width, 720); canvas.height = image.height + 76;
      const context = canvas.getContext('2d'); context.fillStyle = '#ffffff'; context.fillRect(0, 0, canvas.width, canvas.height); context.drawImage(image, (canvas.width - image.width) / 2, 0);
      context.fillStyle = '#233952'; context.font = '14px "Microsoft YaHei", sans-serif';
      context.fillText('InterpreterPLATO · ' + (exportedCaseId ? '实际运行 ' + exportedCaseId + ' · ' : '') + exportedModel.mode, 15, image.height + 27, canvas.width - 30);
      context.fillStyle = '#7d8ca0'; context.font = '11px "Microsoft YaHei", sans-serif';
      context.fillText('当前筛选视图 · 文件关系未独立核验 · PPI 无方向，连通不构成因果证据', 15, image.height + 51, canvas.width - 30);
      canvas.toBlob(blob => { if (blob) download(blob, 'InterpreterPLATO-' + (exportedModel.synthetic ? 'SYNTHETIC-' : exportedModel.payload.run_mode === 'reference_example' ? 'REFERENCE-' : '') + 'network.png'); else error('浏览器未能导出图像。'); }, 'image/png');
    };
    image.onerror = () => error('浏览器未能导出图像。'); image.src = dataUrl;
  }
  function boot() {
    if (!G || !window.InterpreterReference || !window.InterpreterExample || typeof window.cytoscape !== 'function') { error('查看器资源未加载完整。请保留目录内的资源文件，或使用打包后的单文件 HTML。'); return; }
    if (RUNS.length) {
      const buttons = $('run-case-list'); buttons.replaceChildren();
      RUNS.forEach((run, index) => {
        if (run.schema !== 'interpreterplato-local-run/v1') throw new Error('Unknown local run bundle schema');
        const value = index ? 'run:' + index : 'run', p = run.payload, quality = p.network_statistics && p.network_statistics.quality_status;
        $('example-select').insertBefore(new Option('实际运行 · ' + run.source.case_id, value), $('example-select').options[index]);
        const button = element('button', 'run-case-button'); button.type = 'button'; button.dataset.runIndex = String(index); button.setAttribute('aria-pressed', 'false');
        button.append(element('small', '', run.source.case_id), element('strong', '', p.drug1_name + ' + ' + p.drug2_name), element('span', '', (p.cell_line || '未记录细胞系') + ' · 统计质量' + (QUALITY[quality] || quality || '未记录')));
        button.addEventListener('click', () => { $('example-select').value = value; loadBuiltin(); }); buttons.appendChild(button);
      });
      $('run-suite-count').textContent = RUNS.length + ' 个已完成的开发测试实例';
      $('example-select').value = 'run';
    }
    $('show-run-details').addEventListener('click', showRunDetails);
    $('file-input').addEventListener('change', event => openFile(event.target.files[0]));
    $('load-example').addEventListener('click', loadBuiltin); $('example-select').addEventListener('change', loadBuiltin);
    document.querySelectorAll('[data-tour]').forEach(button => button.addEventListener('click', () => tour(button.dataset.tour)));
    $('view-focus').addEventListener('click', () => setReferenceView('focus'));
    $('view-full').addEventListener('click', () => setReferenceView('full'));
    $('show-reference-sources').addEventListener('click', showSources); referenceSources();
    $('result-select').addEventListener('change', event => {
      const index = Number(event.target.value);
      try { error(''); display(index); } catch (exception) {
        const payload = state.rows[index];
        const empty = G.build({ ...payload, evidence_graphs: null });
        empty.payload = payload; empty.warnings = ['证据图未展示：' + exception.message];
        display(index, empty);
        error('该条记录无法展示：' + exception.message);
      }
    });
    document.querySelectorAll('[data-node-group], [data-edge-group]').forEach(input => input.addEventListener('change', applyFilters));
    $('layout-select').addEventListener('change', layout);
    $('toggle-sidebar').addEventListener('click', () => sidebar($('explore-panel').hidden));
    $('toggle-inspector').addEventListener('click', () => inspector($('inspector-panel').hidden));
    $('close-inspector').addEventListener('click', () => { inspector(false); $('toggle-inspector').focus(); });
    $('inspector-panel').addEventListener('keydown', event => { if (event.key === 'Escape') { inspector(false); $('toggle-inspector').focus(); } });
    $('zoom-in').addEventListener('click', () => zoom(1.25));
    $('zoom-out').addEventListener('click', () => zoom(0.8));
    $('fit-graph').addEventListener('click', () => state.cy && state.cy.fit(state.cy.elements().filter(shown), 60));
    $('clear-selection').addEventListener('click', () => { clearHighlight(); $('node-search').value = ''; $('search-status').textContent = '已取消高亮'; });
    $('node-search').addEventListener('input', () => {
      clearHighlight(); const query = $('node-search').value.trim().toLocaleLowerCase();
      if (!query || !state.cy) { $('search-status').textContent = '输入名称可高亮匹配节点'; return; }
      const matches = state.cy.nodes().filter(node => shown(node) && (node.data('label') + ' ' + (node.data('displayLabel') || '') + ' ' + node.data('originalId')).toLocaleLowerCase().includes(query));
      state.cy.elements().addClass('dimmed'); matches.removeClass('dimmed').addClass('highlighted');
      $('search-status').textContent = '匹配 ' + matches.length + ' 个当前可见节点。';
    });
    $('node-search').addEventListener('keydown', event => {
      if (event.key !== 'Enter' || !state.cy) return;
      const matches = state.cy.nodes('.highlighted').filter(shown);
      if (matches.length) { inspect(matches[0]); tab('details'); inspector(true); state.cy.fit(matches, 80); }
    });
    $('export-png').addEventListener('click', () => exportPng().catch(exception => error(exception.message)));
    document.addEventListener('dragover', event => { event.preventDefault(); });
    document.addEventListener('drop', event => { event.preventDefault(); if (event.dataTransfer.files.length !== 1) error('请每次拖入一个解释 JSON 文件。'); else openFile(event.dataTransfer.files[0]); });
    document.querySelectorAll('label.button[for="file-input"]').forEach(label => { label.tabIndex = 0; label.setAttribute('role', 'button'); label.addEventListener('keydown', event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); $('file-input').click(); } }); });
    bindTabs(['details', 'diagnostics'], tab); bindTabs(['paths', 'filters'], exploreTab);
    const resize = () => {
      if (!state.cy) return;
      state.cy.resize();
      if (state.compact !== compactView() && $('layout-select').value === 'preset') layout();
    };
    if (typeof ResizeObserver === 'function') new ResizeObserver(resize).observe($('cy'));
    window.addEventListener('resize', resize);
    loadBuiltin();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot); else boot();
})();
