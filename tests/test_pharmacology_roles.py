"""Drug-specific roles must remain annotations, not new prediction evidence."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json

import networkx as nx
import pytest

from InterpreterPLATO.context import CancerContext, ContextScore
from InterpreterPLATO.mechanism import classify_mechanism, qualified_target_action_status
from InterpreterPLATO.pharmacology import (
    DEFAULT_ROLE_SHA256, PharmacologyRoleCatalog, pharmacology_overview,
    role_composition, summarize_roles, summary_target_label,
)
from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline
from InterpreterPLATO.targets import standardize_smiles

ROOT = Path(__file__).resolve().parents[1]
PRALA = 'C#CCC(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)N[C@@H](CCC(=O)O)C(=O)O'
DOC = 'CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O'


def item(symbol, score=0.5):
    return ContextScore(symbol, symbol, score, 0.35, None, 0.0, 0.0, 0.0,
                        'predicted_only', 'PLATO', {'plato': {'reliable': True}})


def catalog():
    return PharmacologyRoleCatalog(ROOT/'data/pharmacology_roles_v1.json', DEFAULT_ROLE_SHA256)


def write_resource(tmp_path, payload):
    p = tmp_path/'roles.json'
    p.write_text(json.dumps(payload), encoding='utf-8')
    return p, hashlib.sha256(p.read_bytes()).hexdigest()


def resource():
    return json.loads((ROOT/'data/pharmacology_roles_v1.json').read_text(encoding='utf-8'))


def test_bundled_resource_has_portable_lf_bytes_and_pinned_identity():
    raw = (ROOT/'data/pharmacology_roles_v1.json').read_bytes()
    assert b'\r' not in raw, 'The bundled resource must retain LF bytes on every platform'
    assert hashlib.sha256(raw).hexdigest() == DEFAULT_ROLE_SHA256


def test_known_roles_preserve_candidates_scores_order_and_original_evidence():
    original = [item('ATIC', 0.9), item('DHFR', 0.6), item('FPGS', 0.5), item('SLC19A1', 0.4)]
    before = deepcopy(original)
    annotated = catalog().annotate(original, PRALA)
    assert original == before
    for old, new in zip(original, annotated):
        assert old.target_id == new.target_id and old.final_score == new.final_score
        assert old.source_evidence == new.source_evidence
        assert old.evidence_status == new.evidence_status == 'predicted_only'
        assert not old.pharmacology
    groups = summarize_roles(annotated)['groups']
    assert [r['target_id'] for r in groups['direct_pharmacodynamic']] == ['DHFR']
    assert {r['target_id'] for r in groups['drug_handling']} == {'FPGS','SLC19A1'}
    assert [r['target_id'] for r in groups['unresolved']] == ['ATIC']
    assert summary_target_label(annotated) == 'DHFR'
    assert qualified_target_action_status(annotated[1]) == 'unknown'


def test_roles_are_drug_specific_and_family_scope_is_not_isoform_validation():
    original = [item('EGFR', 0.9), item('TUBB3', 0.5), item('ABCB1', 0.4), item('CYP3A4', 0.3)]
    annotated = catalog().annotate(original, DOC)
    summary = summarize_roles(annotated)
    assert not summary['groups']['direct_pharmacodynamic']
    assert [r['target_id'] for r in summary['groups']['family_pharmacodynamic']] == ['TUBB3']
    assert annotated[0].pharmacology['roles'] == ['unknown']
    assert annotated[1].pharmacology['records'][0]['evidence_scope'] == 'target_family'
    assert 'family evidence' in summary_target_label(annotated)
    assert all(t.pharmacology['roles']==['unknown'] for t in catalog().annotate(original, PRALA))
    assert all(t.pharmacology['roles']==['unknown'] for t in catalog().annotate(original, 'CCO'))
    assert not pharmacology_overview(annotated, annotated, ('a','b'))['family_evidence_is_isoform_validation']


def test_structure_stereochemistry_is_not_silently_ignored():
    changed = PRALA.replace('[C@@H]', '[C@H]')
    assert standardize_smiles(changed)['inchikey'] != standardize_smiles(PRALA)['inchikey']
    assert catalog().annotate([item('DHFR')], changed)[0].pharmacology['roles'] == ['unknown']


def test_annotations_never_add_missing_known_target():
    annotated = catalog().annotate([item('ABCB1')], DOC)
    assert [x.target_id for x in annotated] == ['ABCB1']
    assert not summarize_roles(annotated)['groups']['family_pharmacodynamic']
    assert summary_target_label(annotated) == 'no assigned pharmacodynamic target'


def test_binding_to_a_transporter_does_not_make_it_a_drug_handling_role(tmp_path):
    payload = resource()
    r = deepcopy(payload['records'][0])
    r.update(id='SYNTHETIC-SERT',drug_name='Synthetic test compound',drug_inchikey=standardize_smiles('CCO')['inchikey'],target_gene='SLC6A4',process='Synthetic transporter inhibition')
    payload['records'] = [r]
    custom = PharmacologyRoleCatalog(*write_resource(tmp_path,payload))
    annotated = custom.annotate([item('SLC6A4')], 'CCO')
    groups = summarize_roles(annotated)['groups']
    assert len(groups['direct_pharmacodynamic']) == 1 and not groups['drug_handling']


def test_multiple_roles_are_preserved_instead_of_forcing_one_category(tmp_path):
    payload=resource()
    extra=deepcopy(payload['records'][0])
    extra.update(id='SYNTHETIC-MULTIROLE',role='drug_metabolism',action='metabolic_substrate')
    payload['records'].append(extra)
    annotated=PharmacologyRoleCatalog(*write_resource(tmp_path,payload)).annotate([item('DHFR')],PRALA)
    groups=summarize_roles(annotated)['groups']
    assert len(groups['direct_pharmacodynamic'])==len(groups['drug_handling'])==1
    assert role_composition(['DHFR'],annotated)['drug_handling_only'] is False


def test_pathway_counts_do_not_promote_adme_or_unknown_to_pharmacodynamics():
    profile=catalog().annotate([item('ABCB1'),item('CYP3A4'),item('SLCO1B3'),item('EGFR')],DOC)
    assert role_composition(['ABCB1','CYP3A4','SLCO1B3'],profile)['drug_handling_only'] is True
    for genes in ([],['EGFR'],['ABCB1','EGFR'],['NOT_IN_PROFILE']):
        r=role_composition(genes,profile)
        assert r['drug_handling_only'] is False
        assert r['pharmacodynamic_enrichment_tested'] is False


@pytest.mark.parametrize('mutation', ['missing_ref','bad_role','association_as_target','duplicate','wildcard','partial_inchikey'])
def test_invalid_role_resources_fail_closed(tmp_path,mutation):
    payload=resource()
    r=payload['records'][0]
    if mutation=='missing_ref':r['reference_ids']=['NO_REFERENCE']
    elif mutation=='bad_role':r['role']='target_because_highly_expressed'
    elif mutation=='association_as_target':r['evidence_scope']='association_only'
    elif mutation=='duplicate':payload['records'].append(deepcopy(r))
    elif mutation=='wildcard':r['target_gene']='TUBB*'
    elif mutation=='partial_inchikey':r['drug_inchikey']=r['drug_inchikey'].split('-')[0]
    with pytest.raises(ValueError):PharmacologyRoleCatalog(*write_resource(tmp_path,payload))


def test_corrupted_resource_and_missing_hash_are_rejected(tmp_path):
    path,sha=write_resource(tmp_path,resource())
    with pytest.raises(ValueError):PharmacologyRoleCatalog(path,'')
    path.write_text(path.read_text()+' ',encoding='utf-8')
    with pytest.raises(ValueError):PharmacologyRoleCatalog(path,sha)


def test_role_annotations_leave_core_classification_and_statistics_identical():
    context=CancerContext()
    symbols=['DHFR','FPGS','SLC19A1','TUBB3','ABCB1','CYP3A4','EGFR']
    graph=nx.cycle_graph(symbols+[f'G{i}' for i in range(80)])
    nx.set_edge_attributes(graph,1.0,'distance')
    context.load_ppi(graph);context.default_disease_genes={'DHFR'}
    left=[item('DHFR'),item('FPGS'),item('SLC19A1')]
    right=[item('TUBB3'),item('ABCB1'),item('CYP3A4'),item('EGFR')]
    before=classify_mechanism(left,right,context,network_randomizations=1000)
    after=classify_mechanism(catalog().annotate(left,PRALA),catalog().annotate(right,DOC),context,network_randomizations=1000)
    assert before.to_dict()==after.to_dict()


def test_pipeline_role_toggle_changes_narrative_only_and_audits_resource(tmp_path):
    cache=tmp_path/'cache';cache.mkdir()
    def rows(symbols):
        return [{'target_symbol':s,'target_name':s,'organism':'Homo sapiens','taxid':'9606','target_type':'single protein','score':12.0,'reliable':'yes','max_similarity':0.8} for s in symbols]
    hashes=[archive_verified_plato_cache(cache,PRALA,'a',rows(['DHFR','FPGS','SLC19A1'])),archive_verified_plato_cache(cache,DOC,'b',rows(['TUBB3','ABCB1','EGFR']))]
    pipe=strict_pipeline(cache,hashes)
    pipe.load_demo_resources()
    pipe.config['network']['randomizations']=20
    pipe.config['credibility']['stability_randomizations']=20
    enabled=pipe.explain_combination('name-does-not-control-matching','another-alias',PRALA,DOC,advanced=True).to_dict()
    pipe.config['pharmacology_roles']['enabled']=False
    disabled=pipe.explain_combination('name-does-not-control-matching','another-alias',PRALA,DOC,advanced=True).to_dict()
    assert enabled['network_statistics']==disabled['network_statistics']
    for k in ['mechanism_confidence','mechanism_type','evidence_credibility_score','validation']:
        assert enabled[k]==disabled[k]
    assert enabled['mechanism']['description'] != disabled['mechanism']['description']
    assert enabled['audit']['pharmacology_roles']['resource_sha256']==DEFAULT_ROLE_SHA256
    assert disabled['audit']['pharmacology_roles']['enabled'] is False
    assert enabled['pharmacology']['drug1']['groups']['direct_pharmacodynamic'][0]['target_id']=='DHFR'
    assert enabled['pharmacology']['drug2']['groups']['family_pharmacodynamic'][0]['target_id']=='TUBB3'
    assert 'family evidence' in enabled['integrated_explanation']['summary']
    assert enabled['drug1_targets']['top_targets_semantics']=='context_rank_not_primary_molecular_target_rank'
