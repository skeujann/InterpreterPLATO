"""Path correctness, evidence-graph integrity, and API/CLI wiring."""

import json
import random
import shutil
import sys

import networkx as nx
import pytest

from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.advanced import build_evidence_graph, local_evidence_paths
from InterpreterPLATO.connectivity import bounded_weighted_paths, cross_drug_paths, validate_path_search
from InterpreterPLATO.context import CancerContext, ContextScore, edge_distance
from InterpreterPLATO.scripts import query_pair_and_explain, run_interpreter
from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline


def target(symbol, score=0.7):
    return ContextScore(symbol, symbol, score, score, None, 0.0, 0.0, 0.0,
                        "predicted_only", "PLATO", {})


def context(edges=(), disease=("D",)):
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges, weight="distance")
    result = CancerContext()
    result.load_ppi(graph)
    result.default_disease_genes = set(disease)
    return result


def test_bounded_search_recovers_more_expensive_short_route():
    ctx = context([("S", "A", 1), ("A", "B", 1), ("B", "D", 1), ("S", "D", 7)])
    assert nx.shortest_path(ctx.ppi, "S", "D", weight="distance") == ["S", "A", "B", "D"]
    assert bounded_weighted_paths(ctx.ppi, ["S"], ["D"], 2) == {"S": ["S", "D"]}


def test_bounded_search_keeps_shorter_prefix_even_when_it_costs_more():
    ctx = context([("D", "A", 1), ("A", "B", 1), ("B", "X", 1), ("D", "X", 6), ("X", "S", 1)])
    assert bounded_weighted_paths(ctx.ppi, ["S"], ["D"], 3)["S"] == ["S", "X", "D"]


@pytest.mark.parametrize("seed", range(12))
def test_bounded_search_agrees_with_exhaustive_simple_paths(seed):
    rng = random.Random(seed)
    graph = nx.Graph()
    nodes = list("ABCDEFG")
    graph.add_nodes_from(nodes)
    for index, left in enumerate(nodes):
        for right in nodes[index + 1:]:
            if rng.random() < 0.4:
                graph.add_edge(left, right, distance=rng.randint(1, 5))
    destinations = ["F", "G"]
    found = bounded_weighted_paths(graph, nodes, destinations, 4)
    for source in nodes:
        routes = [[source]] if source in destinations else []
        for destination in destinations:
            if destination != source:
                routes.extend(nx.all_simple_paths(graph, source, destination, cutoff=4))
        if not routes:
            assert source not in found
        else:
            expected = min(routes, key=lambda route: (
                sum(edge_distance(a, b, graph[a][b]) for a, b in zip(route, route[1:])),
                len(route), tuple(route),
            ))
            assert found[source] == expected


def test_tied_routes_do_not_depend_on_edge_or_destination_insertion_order():
    edges = [("S", "B", 1), ("B", "D", 1), ("S", "A", 1), ("A", "E", 1)]
    one, two = context(edges), context(list(reversed(edges)))
    assert bounded_weighted_paths(one.ppi, ["S"], ["D", "E"], 2) == {"S": ["S", "A", "E"]}
    assert bounded_weighted_paths(one.ppi, ["S"], ["D", "E"], 2) == bounded_weighted_paths(two.ppi, ["S"], ["E", "D"], 2)


def test_diagnostics_distinguish_mapping_disconnection_length_and_display_limits():
    ctx = context([("A", "D", 1), ("B", "D", 1), ("F", "Y", 1), ("Y", "Z", 1), ("Z", "D", 1), ("X", "Q", 1)])
    diagnostics = {}
    paths = local_evidence_paths("a", "b", [target("MISSING"), target("X"), target("F"), target("A"), target("B", .6)], [], ctx, top_k=1, max_length=2, diagnostics=diagnostics)
    assert len(paths) == 1 and paths[0]["source_target"] == "A"
    statuses = {item["target"]: item["status"] for item in diagnostics["drugs"][0]["targets"]}
    assert statuses == {"MISSING": "target_not_in_ppi", "X": "disconnected_from_disease", "F": "path_length_limit", "A": "shown", "B": "display_limit"}
    assert diagnostics["drugs"][0]["eligible_target_count"] == 5
    assert diagnostics["candidate_path_count"] == 2


def test_both_drugs_have_display_slots_and_zero_budget_does_not_mean_disconnected():
    ctx = context([("A", "D", 1), ("B", "D", 1), ("C", "D", 1)])
    profiles = ([target("A", .9), target("B", .8)], [target("C", .1)])
    paths = local_evidence_paths("a", "b", *profiles, ctx, top_k=2)
    assert {item["drug"] for item in paths} == {"a", "b"}
    diagnostics = {}
    assert local_evidence_paths("a", "b", *profiles, ctx, top_k=0, diagnostics=diagnostics) == []
    assert diagnostics["candidate_path_count"] == 3
    assert all(item["status"] == "display_limit" for drug in diagnostics["drugs"] for item in drug["targets"])


def test_cross_drug_paths_work_without_disease_and_preserve_graph():
    ctx = context([("A", "BRIDGE", 1), ("BRIDGE", "B", 1)], disease=())
    before = nx.node_link_data(ctx.ppi)
    assert local_evidence_paths("a", "b", [target("A")], [target("B")], ctx) == []
    paths, report = cross_drug_paths(("a", "b"), ([target("A")], [target("B")]), ctx, max_paths=5, max_length=2)
    assert report["status"] == "connected"
    assert [item["path"] for item in paths] == [["A", "BRIDGE", "B"]]
    assert nx.node_link_data(ctx.ppi) == before
    graph = build_evidence_graph("a", "b", [target("A")], [target("B")], ctx, paths, [])
    ppi_edges = [edge for edge in graph["edges"] if edge["type"] == "PPI"]
    assert len(ppi_edges) == 2
    assert all(edge["directed"] is False for edge in ppi_edges)


def test_zero_hop_is_membership_not_ppi_support_and_bridge_disease_nodes_are_linked():
    ctx = context([("A", "D", 1)], disease=("D",))
    paths = local_evidence_paths("a", "b", [target("A")], [target("D")], ctx)
    direct = next(item for item in paths if item["source_target"] == "D")
    assert direct["direct_disease_membership"] is True
    assert direct["edge_confidence_coverage"] == 0.0
    assert direct["ppi_edge_count"] == 0
    assert direct["edge_confidence_applicable"] is False
    # D is introduced solely as the path endpoint in this graph call.
    graph = build_evidence_graph("a", "b", [target("A")], [], ctx, paths, [])
    assert {"source": "gene:D", "target": "disease_module:selected_context", "type": "belongs_to_disease_module"} in graph["edges"]
    node_ids = {item["id"] for item in graph["nodes"]}
    assert all(edge["source"] in node_ids and edge["target"] in node_ids for edge in graph["edges"])


def test_only_selected_disease_module_is_attached_and_ppi_edges_are_deduplicated():
    ctx = context([("A", "D", 1)])
    ctx.disease_sets = {"chosen": {"D"}, "unrelated": {"A"}}
    paths = [{"path": ["A", "D"], "score": .5}, {"path": ["D", "A"], "score": .8, "path_kind": "cross_drug_target"}]
    graph = build_evidence_graph("a", "b", [target("A")], [], ctx, paths, [], disease_label="chosen")
    assert [node["id"] for node in graph["nodes"] if node["type"] == "disease_module"] == ["disease_module:chosen"]
    edges = [edge for edge in graph["edges"] if edge["type"] == "PPI"]
    assert len(edges) == 1 and edges[0]["path_score"] == .8
    assert edges[0]["path_kinds"] == ["cross_drug_target", "target_to_disease"]


def test_evidence_graph_rejects_invented_edge():
    ctx = context([("A", "X", 1), ("D", "Y", 1)])
    with pytest.raises(ValueError, match="not in the loaded PPI"):
        build_evidence_graph("a", "b", [target("A")], [target("D")], ctx, [{"path": ["A", "D"]}], [])


@pytest.mark.parametrize("edges,disease,status", [
    ([], ["D"], "missing_ppi"),
    ([("A", "B", 1)], [], "missing_disease_context"),
    ([("A", "B", 1)], ["D"], "disease_not_in_ppi"),
])
def test_missing_resources_are_reported_without_inventing_routes(edges, disease, status):
    ctx = context(edges, disease)
    diagnostics = {}
    assert local_evidence_paths("a", "b", [target("A")], [], ctx, diagnostics=diagnostics) == []
    assert diagnostics["drugs"][0]["targets"][0]["status"] == status


@pytest.mark.parametrize("left,right,max_length,status", [
    (["A"], ["D"], 5, "disconnected_target_sets"),
    (["A"], ["MISSING"], 5, "target_mapping_incomplete"),
    (["A"], ["B"], 0, "path_length_limit"),
])
def test_cross_drug_negative_cases_do_not_create_connections(left, right, max_length, status):
    ctx = context([("A", "B", 1), ("D", "E", 1)])
    paths, report = cross_drug_paths(("a", "b"), ([target(s) for s in left], [target(s) for s in right]), ctx, max_paths=5, max_length=max_length)
    assert paths == []
    assert report["status"] == status


def test_shared_target_is_zero_hop_and_does_not_invent_a_self_edge():
    ctx = context([("A", "B", 1)])
    paths, report = cross_drug_paths(("a", "b"), ([target("A")], [target("A")]), ctx, max_paths=5, max_length=0)
    assert paths[0]["shared_target"] is True and paths[0]["path"] == ["A"]
    assert report["shared_targets"] == ["A"]
    graph = build_evidence_graph("a", "b", [target("A")], [target("A")], ctx, paths, [])
    assert all(edge["type"] != "PPI" for edge in graph["edges"])


@pytest.mark.parametrize("options", [{"max_length": -1}, {"max_length": True}, {"max_paths": 1.2}, {"max_paths": 101}, {"unknown": 2}, None])
def test_configuration_rejects_invalid_values(options):
    with pytest.raises(ValueError):
        validate_path_search(options)
    with pytest.raises(ValueError):
        ExplanationPipeline(config={"path_search": options})


def cached_pair(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    def rows(symbol):
        return [{"target_symbol": symbol, "target_name": symbol, "organism": "Homo sapiens", "taxid": "9606", "target_type": "single protein", "score": 12., "reliable": "yes", "max_similarity": .8}]
    hashes = [archive_verified_plato_cache(cache, "CCO", "a", rows("EGFR")), archive_verified_plato_cache(cache, "CCN", "b", rows("MTOR"))]
    return cache, hashes


def test_pipeline_wires_options_serialization_and_disabled_status_without_statistical_bonus(tmp_path):
    cache, hashes = cached_pair(tmp_path)
    pipe = strict_pipeline(cache, hashes, advanced=True, path_search={"max_length": 0, "max_paths": 1, "max_cross_drug_paths": 0})
    pipe.load_demo_resources()
    # This regression must not rely on unversioned local demo data.
    pipe.context.load_ppi(nx.Graph([("EGFR", "MTOR")]))
    restricted = pipe.explain_combination("a", "b", "CCO", "CCN").to_dict()
    assert restricted["network_connectivity"]["max_length"] == 0
    assert len(restricted["top_evidence_paths"]) <= 1
    assert restricted["local_evidence"]["cross_drug_paths"] == []
    pipe.config["path_search"] = {"max_length": 5, "max_paths": 10, "max_cross_drug_paths": 5}
    expanded = pipe.explain_combination("a", "b", "CCO", "CCN").to_dict()
    assert expanded["local_evidence"]["cross_drug_paths"]
    assert expanded["network_connectivity"] == expanded["local_evidence"]["connectivity"]
    assert restricted["network_statistics"] == expanded["network_statistics"]
    assert restricted["mechanism_type"] == expanded["mechanism_type"]
    assert restricted["mechanism_confidence"] == expanded["mechanism_confidence"]
    assert restricted["evidence_credibility_score"] == expanded["evidence_credibility_score"]
    output = tmp_path / "explanation.json"
    pipe.save_result(pipe.explain_combination("a", "b", "CCO", "CCN", advanced=False), output)
    assert json.loads(output.read_text(encoding="utf-8"))["network_connectivity"] == {"status": "disabled"}
    assert (tmp_path / "explanation.json.sha256").is_file()
    pipe.config["path_search"]["max_length"] = -1
    with pytest.raises(ValueError, match="path_search.max_length"):
        pipe.explain_combination("a", "b", "CCO", "CCN")


@pytest.mark.parametrize("command", ["summary", "batch"])
def test_cli_options_reach_saved_report(tmp_path, monkeypatch, command):
    cache, hashes = cached_pair(tmp_path)
    output = tmp_path / "cli.json"
    argv = ["run_interpreter", command, "--demo", "--quick", "--plato-cache", str(cache), "--output", str(output), "--path-max-length", "0", "--max-evidence-paths", "1", "--max-cross-drug-paths", "0"]
    for digest in hashes:
        argv += ["--trusted-plato-sha256", digest]
    if command == "summary":
        argv += ["--drug1", "a", "--drug2", "b", "--drug1-smiles", "CCO", "--drug2-smiles", "CCN"]
    else:
        batch = tmp_path / "pairs.json"
        batch.write_text(json.dumps([{"drug1_name": "a", "drug2_name": "b", "drug1_smiles": "CCO", "drug2_smiles": "CCN"}]), encoding="utf-8")
        argv += ["--input", str(batch)]
    # Keep the real CLI flow, but mark these synthetic test inputs explicitly.
    def fixture_pipeline(config):
        return ExplanationPipeline(config={**config, "run_mode": "test", "generated_from_test_fixture": True})
    monkeypatch.setattr(run_interpreter, "ExplanationPipeline", fixture_pipeline)
    monkeypatch.setattr(sys, "argv", argv)
    assert run_interpreter.main() == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    if command == "batch":
        payload = payload[0]
    assert payload["network_connectivity"]["max_length"] == 0
    assert payload["network_connectivity"]["max_paths"] == 1
    assert payload["network_connectivity"]["cross_drug"]["max_paths"] == 0
    assert payload["local_evidence"]["cross_drug_paths"] == []


def test_query_pair_options_reach_manifest_and_explanation_without_network(tmp_path, monkeypatch):
    cache, hashes = cached_pair(tmp_path)
    run_dir = tmp_path / "pair_run"
    called = []
    def fake_capture(**options):
        called.append(options["drug_name"])
        shutil.copytree(cache, options["cache_dir"], dirs_exist_ok=True)
        return {"raw_sha256": hashes[0 if options["drug_name"] == "a" else 1], "synthetic_fixture": True}
    def fixture_pipeline(config):
        return ExplanationPipeline(config={**config, "run_mode": "test", "generated_from_test_fixture": True})
    argv = ["query_pair", "--drug1", "a", "--drug2", "b", "--drug1-smiles", "CCO", "--drug2-smiles", "CCN", "--demo", "--quick", "--delay-seconds", "0", "--run-dir", str(run_dir), "--path-max-length", "0", "--max-evidence-paths", "1", "--max-cross-drug-paths", "0"]
    for digest in hashes:
        argv += ["--trusted-plato-sha256", digest]
    monkeypatch.setattr(query_pair_and_explain, "capture_plato_single_query", fake_capture)
    monkeypatch.setattr(query_pair_and_explain, "ExplanationPipeline", fixture_pipeline)
    monkeypatch.setattr(sys, "argv", argv)
    assert query_pair_and_explain.main() == 0
    assert called == ["a", "b"]
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    payload = json.loads((run_dir / "explanation.json").read_text(encoding="utf-8"))
    assert manifest["path_search"] == {"max_length": 0, "max_paths": 1, "max_cross_drug_paths": 0}
    assert payload["network_connectivity"]["max_length"] == 0
    assert payload["local_evidence"]["cross_drug_paths"] == []


def test_query_pair_rejects_invalid_options_before_external_requests(monkeypatch):
    def forbidden_capture(**options):
        pytest.fail("Invalid path options must fail before querying PLATO")
    monkeypatch.setattr(query_pair_and_explain, "capture_plato_single_query", forbidden_capture)
    monkeypatch.setattr(sys, "argv", ["query_pair", "--drug1", "a", "--drug2", "b", "--path-max-length", "-1"])
    with pytest.raises(SystemExit) as error:
        query_pair_and_explain.main()
    assert error.value.code == 2
