/* Public reference annotations. See REFERENCE_CASE.md for attribution and scope. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.InterpreterReference = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  return {
  "result_schema_version": "5.4",
  "drug1_name": "Dabrafenib",
  "drug2_name": "Trametinib",
  "run_mode": "reference_example",
  "generated_from_test_fixture": false,
  "formal_result_eligible": false,
  "predicted_synergy": null,
  "mechanism_type": "not_computed",
  "metadata": {
    "reference_example": true,
    "calculation_performed": false,
    "prepared_on": "2026-10-21",
    "title": "BRAF / MEK 双药文献导览",
    "context": "BRAF V600 突变黑色素瘤 · 文献背景，非细胞系实验",
    "notice": "人工整理的公共注释实例，非 PLATO / InterpreterPLATO 计算结果。无协同分数、无富集检验、无患者或细胞系推断。",
    "selection_policy": "Seven preselected canonical human UniProt accessions related to MAPK; retain all unique non-self, same-species INTERACTION pairs with exact accession matches. Isoforms are not collapsed into canonical accessions. Paths are editorially selected existing-edge routes, not an exhaustive search.",
    "visualization_positions": {
      "drug:Dabrafenib": [
        -320,
        -150
      ],
      "drug:Trametinib": [
        -320,
        40
      ],
      "gene:MAP2K1": [
        0,
        -30
      ],
      "gene:MAP2K2": [
        0,
        110
      ],
      "gene:BRAF": [
        320,
        -50
      ],
      "gene:RAF1": [
        100,
        -230
      ],
      "gene:KSR1": [
        -120,
        260
      ],
      "gene:MAPK1": [
        130,
        260
      ],
      "gene:MAPK3": [
        330,
        260
      ],
      "pathway:R-HSA-5673001": [
        560,
        240
      ],
      "disease_module:literature-context": [
        665,
        -50
      ]
    },
    "references": [
      {
        "title": "组合研究 · Flaherty et al., NEJM 2012",
        "url": "https://pubmed.ncbi.nlm.nih.gov/23020132/",
        "detail": "PMID 23020132 · DOI 10.1056/NEJMoa1210093。只引用组合与研究背景，不转载正文。"
      },
      {
        "title": "药物靶点 · NCI",
        "url": "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/dabrafenib",
        "detail": "Dabrafenib 的 BRAF 药理靶点注释。"
      },
      {
        "title": "MEK1 / MEK2 · NCI",
        "url": "https://www.cancer.gov/publications/dictionaries/cancer-terms/def/mekinist",
        "detail": "Trametinib 的 MEK1 / MEK2 药理靶点注释。"
      },
      {
        "title": "人类蛋白互作 · UniProtKB / IntAct",
        "url": "https://www.uniprot.org/uniprotkb/P15056/entry",
        "detail": "7 份 UniProt 条目的 INTERACTION 子集，逐边记录 IntAct 标识、原始响应 SHA-256 与条目版本。实验计数不是置信度。"
      },
      {
        "title": "通路成员 · Reactome",
        "url": "https://reactome.org/content/detail/R-HSA-5673001",
        "detail": "R-HSA-5673001 的 7 个选定参与蛋白；不是通路富集结果。"
      }
    ],
    "licenses": [
      {
        "resource": "UniProt Consortium — selected UniProtKB annotations",
        "license": "CC BY 4.0",
        "url": "https://creativecommons.org/licenses/by/4.0/",
        "source": "https://www.uniprot.org/help/license",
        "modifications": "Selected identifiers and interaction facts; deduplicated undirected pairs; translated display labels. No source prose reproduced."
      },
      {
        "resource": "Reactome — selected pathway participant facts",
        "license": "CC0",
        "url": "https://creativecommons.org/publicdomain/zero/1.0/",
        "source": "https://reactome.org/license",
        "modifications": "Selected seven participant identifiers; not a database dump or copied diagram."
      }
    ],
    "focus_view": {
      "endpoint": "disease_module:literature-context",
      "nodes": [
        "drug:Dabrafenib",
        "drug:Trametinib",
        "gene:MAP2K1",
        "gene:MAP2K2",
        "gene:BRAF",
        "disease_module:literature-context"
      ],
      "semantics": "Reading order only: drugs, targets and PPI, BRAF, disease-module literature context. The terminal node is a single-gene literature annotation, not a computed disease gene set. Position and emphasis do not imply shared direct targeting, biological direction or causality.",
      "anchor": "gene:BRAF"
    }
  },
  "drug1_targets": {
    "targets": [
      {
        "target_id": "BRAF",
        "evidence_status": "literature_annotation"
      }
    ]
  },
  "drug2_targets": {
    "targets": [
      {
        "target_id": "MAP2K1",
        "evidence_status": "literature_annotation"
      },
      {
        "target_id": "MAP2K2",
        "evidence_status": "literature_annotation"
      }
    ]
  },
  "evidence_graphs": {
    "explanation_evidence_graph": {
      "nodes": [
        {
          "id": "drug:Dabrafenib",
          "type": "drug",
          "name": "Dabrafenib",
          "source_url": "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/dabrafenib"
        },
        {
          "id": "drug:Trametinib",
          "type": "drug",
          "name": "Trametinib",
          "source_url": "https://www.cancer.gov/publications/dictionaries/cancer-terms/def/mekinist"
        },
        {
          "id": "gene:BRAF",
          "type": "gene",
          "name": "BRAF",
          "alias": "BRAF",
          "uniprot_accession": "P15056",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/P15056/entry"
        },
        {
          "id": "gene:MAP2K1",
          "type": "gene",
          "name": "MAP2K1 · MEK1",
          "alias": "MEK1",
          "uniprot_accession": "Q02750",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/Q02750/entry"
        },
        {
          "id": "gene:MAP2K2",
          "type": "gene",
          "name": "MAP2K2 · MEK2",
          "alias": "MEK2",
          "uniprot_accession": "P36507",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/P36507/entry"
        },
        {
          "id": "gene:MAPK1",
          "type": "gene",
          "name": "MAPK1 · ERK2",
          "alias": "ERK2",
          "uniprot_accession": "P28482",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/P28482/entry"
        },
        {
          "id": "gene:MAPK3",
          "type": "gene",
          "name": "MAPK3 · ERK1",
          "alias": "ERK1",
          "uniprot_accession": "P27361",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/P27361/entry"
        },
        {
          "id": "gene:RAF1",
          "type": "gene",
          "name": "RAF1",
          "alias": "RAF1",
          "uniprot_accession": "P04049",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/P04049/entry"
        },
        {
          "id": "gene:KSR1",
          "type": "gene",
          "name": "KSR1",
          "alias": "KSR1",
          "uniprot_accession": "Q8IVT5",
          "taxon_id": 9606,
          "source_url": "https://www.uniprot.org/uniprotkb/Q8IVT5/entry"
        },
        {
          "id": "pathway:R-HSA-5673001",
          "type": "pathway",
          "name": "RAF/MAP kinase",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "membership_semantics": "Selected Reactome reference-entity participants; not enrichment or a complete pathway."
        },
        {
          "id": "disease_module:literature-context",
          "type": "disease_module",
          "name": "BRAF V600 melanoma",
          "source_url": "https://pubmed.ncbi.nlm.nih.gov/23020132/",
          "membership_semantics": "Manually curated trial-context annotation of BRAF only; not a COSMIC or computed disease gene set."
        }
      ],
      "edges": [
        {
          "source": "gene:BRAF",
          "target": "gene:MAP2K1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P15056",
              "entry_version": 269,
              "source_url": "https://rest.uniprot.org/uniprotkb/P15056.json",
              "source_sha256": "f4d1e9e4534be377e91691c09016add1972ecb6f149195320fad2b40260cdc46",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P15056",
                  "intActId": "EBI-365980"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q02750",
                  "geneName": "MAP2K1",
                  "intActId": "EBI-492564"
                },
                "numberOfExperiments": 65,
                "organismDiffer": false
              }
            },
            {
              "entry": "Q02750",
              "entry_version": 249,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q02750.json",
              "source_sha256": "cf5f54050e2a6f837987900ff7145b48bf4ca7ab1a4e49d2710889bd9e59ba24",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q02750",
                  "intActId": "EBI-492564"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P15056",
                  "geneName": "BRAF",
                  "intActId": "EBI-365980"
                },
                "numberOfExperiments": 65,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:BRAF",
          "target": "gene:MAP2K2",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P15056",
              "entry_version": 269,
              "source_url": "https://rest.uniprot.org/uniprotkb/P15056.json",
              "source_sha256": "f4d1e9e4534be377e91691c09016add1972ecb6f149195320fad2b40260cdc46",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P15056",
                  "intActId": "EBI-365980"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P36507",
                  "geneName": "MAP2K2",
                  "intActId": "EBI-1056930"
                },
                "numberOfExperiments": 13,
                "organismDiffer": false
              }
            },
            {
              "entry": "P36507",
              "entry_version": 254,
              "source_url": "https://rest.uniprot.org/uniprotkb/P36507.json",
              "source_sha256": "d316777cec1ef826a0ed062f7fc521d9729c782df917a46066331e00949186de",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P36507",
                  "intActId": "EBI-1056930"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P15056",
                  "geneName": "BRAF",
                  "intActId": "EBI-365980"
                },
                "numberOfExperiments": 13,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:BRAF",
          "target": "gene:RAF1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P15056",
              "entry_version": 269,
              "source_url": "https://rest.uniprot.org/uniprotkb/P15056.json",
              "source_sha256": "f4d1e9e4534be377e91691c09016add1972ecb6f149195320fad2b40260cdc46",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P15056",
                  "intActId": "EBI-365980"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P04049",
                  "geneName": "RAF1",
                  "intActId": "EBI-365996"
                },
                "numberOfExperiments": 73,
                "organismDiffer": false
              }
            },
            {
              "entry": "P04049",
              "entry_version": 271,
              "source_url": "https://rest.uniprot.org/uniprotkb/P04049.json",
              "source_sha256": "6e8cc60c356b3b0578a3c04946af28ecd31aee42960929870baab327dc6697ca",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P04049",
                  "intActId": "EBI-365996"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P15056",
                  "geneName": "BRAF",
                  "intActId": "EBI-365980"
                },
                "numberOfExperiments": 73,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:KSR1",
          "target": "gene:MAP2K2",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P36507",
              "entry_version": 254,
              "source_url": "https://rest.uniprot.org/uniprotkb/P36507.json",
              "source_sha256": "d316777cec1ef826a0ed062f7fc521d9729c782df917a46066331e00949186de",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P36507",
                  "intActId": "EBI-1056930"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q8IVT5",
                  "geneName": "KSR1",
                  "intActId": "EBI-486984"
                },
                "numberOfExperiments": 10,
                "organismDiffer": false
              }
            },
            {
              "entry": "Q8IVT5",
              "entry_version": 197,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q8IVT5.json",
              "source_sha256": "d824c8c065adaa9f2521522d137796d779b5e9cc16ce153948f2142ae9b8b22e",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q8IVT5",
                  "intActId": "EBI-486984"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P36507",
                  "geneName": "MAP2K2",
                  "intActId": "EBI-1056930"
                },
                "numberOfExperiments": 10,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:KSR1",
          "target": "gene:RAF1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P04049",
              "entry_version": 271,
              "source_url": "https://rest.uniprot.org/uniprotkb/P04049.json",
              "source_sha256": "6e8cc60c356b3b0578a3c04946af28ecd31aee42960929870baab327dc6697ca",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P04049",
                  "intActId": "EBI-365996"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q8IVT5",
                  "geneName": "KSR1",
                  "intActId": "EBI-486984"
                },
                "numberOfExperiments": 5,
                "organismDiffer": false
              }
            },
            {
              "entry": "Q8IVT5",
              "entry_version": 197,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q8IVT5.json",
              "source_sha256": "d824c8c065adaa9f2521522d137796d779b5e9cc16ce153948f2142ae9b8b22e",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q8IVT5",
                  "intActId": "EBI-486984"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P04049",
                  "geneName": "RAF1",
                  "intActId": "EBI-365996"
                },
                "numberOfExperiments": 5,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:MAP2K1",
          "target": "gene:MAPK1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "Q02750",
              "entry_version": 249,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q02750.json",
              "source_sha256": "cf5f54050e2a6f837987900ff7145b48bf4ca7ab1a4e49d2710889bd9e59ba24",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q02750",
                  "intActId": "EBI-492564"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P28482",
                  "geneName": "MAPK1",
                  "intActId": "EBI-959949"
                },
                "numberOfExperiments": 2,
                "organismDiffer": false
              }
            },
            {
              "entry": "P28482",
              "entry_version": 261,
              "source_url": "https://rest.uniprot.org/uniprotkb/P28482.json",
              "source_sha256": "5166e9247fd19dd7b92dde5db62e68d67d0db808ee9582c078d33ca64ac3a54b",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P28482",
                  "intActId": "EBI-959949"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q02750",
                  "geneName": "MAP2K1",
                  "intActId": "EBI-492564"
                },
                "numberOfExperiments": 2,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:MAP2K1",
          "target": "gene:MAPK3",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "Q02750",
              "entry_version": 249,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q02750.json",
              "source_sha256": "cf5f54050e2a6f837987900ff7145b48bf4ca7ab1a4e49d2710889bd9e59ba24",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q02750",
                  "intActId": "EBI-492564"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P27361",
                  "geneName": "MAPK3",
                  "intActId": "EBI-73995"
                },
                "numberOfExperiments": 2,
                "organismDiffer": false
              }
            },
            {
              "entry": "P27361",
              "entry_version": 258,
              "source_url": "https://rest.uniprot.org/uniprotkb/P27361.json",
              "source_sha256": "34ab28a7c926756d4461e5486239b127c95382ef52515b927296658386caea1a",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P27361",
                  "intActId": "EBI-73995"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q02750",
                  "geneName": "MAP2K1",
                  "intActId": "EBI-492564"
                },
                "numberOfExperiments": 2,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:MAP2K1",
          "target": "gene:RAF1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "Q02750",
              "entry_version": 249,
              "source_url": "https://rest.uniprot.org/uniprotkb/Q02750.json",
              "source_sha256": "cf5f54050e2a6f837987900ff7145b48bf4ca7ab1a4e49d2710889bd9e59ba24",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "Q02750",
                  "intActId": "EBI-492564"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P04049",
                  "geneName": "RAF1",
                  "intActId": "EBI-365996"
                },
                "numberOfExperiments": 39,
                "organismDiffer": false
              }
            },
            {
              "entry": "P04049",
              "entry_version": 271,
              "source_url": "https://rest.uniprot.org/uniprotkb/P04049.json",
              "source_sha256": "6e8cc60c356b3b0578a3c04946af28ecd31aee42960929870baab327dc6697ca",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P04049",
                  "intActId": "EBI-365996"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "Q02750",
                  "geneName": "MAP2K1",
                  "intActId": "EBI-492564"
                },
                "numberOfExperiments": 39,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "gene:MAP2K2",
          "target": "gene:RAF1",
          "type": "PPI",
          "directed": false,
          "source_dataset": "UniProtKB INTERACTION (IntAct-linked annotation)",
          "confidence": null,
          "confidence_available": false,
          "context_limit": "Database interaction annotation; not necessarily direct binding, variant-specific, melanoma-specific, or evidence for drug synergy.",
          "evidence": [
            {
              "entry": "P36507",
              "entry_version": 254,
              "source_url": "https://rest.uniprot.org/uniprotkb/P36507.json",
              "source_sha256": "d316777cec1ef826a0ed062f7fc521d9729c782df917a46066331e00949186de",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P36507",
                  "intActId": "EBI-1056930"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P04049",
                  "geneName": "RAF1",
                  "intActId": "EBI-365996"
                },
                "numberOfExperiments": 7,
                "organismDiffer": false
              }
            },
            {
              "entry": "P04049",
              "entry_version": 271,
              "source_url": "https://rest.uniprot.org/uniprotkb/P04049.json",
              "source_sha256": "6e8cc60c356b3b0578a3c04946af28ecd31aee42960929870baab327dc6697ca",
              "interaction_record": {
                "interactantOne": {
                  "uniProtKBAccession": "P04049",
                  "intActId": "EBI-365996"
                },
                "interactantTwo": {
                  "uniProtKBAccession": "P36507",
                  "geneName": "MAP2K2",
                  "intActId": "EBI-1056930"
                },
                "numberOfExperiments": 7,
                "organismDiffer": false
              }
            }
          ]
        },
        {
          "source": "drug:Dabrafenib",
          "target": "gene:BRAF",
          "type": "known_target",
          "evidence_status": "literature_annotation",
          "source_url": "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/dabrafenib",
          "scope": "Selected documented pharmacological targets; not an exhaustive target profile or PLATO output."
        },
        {
          "source": "drug:Trametinib",
          "target": "gene:MAP2K1",
          "type": "known_target",
          "evidence_status": "literature_annotation",
          "source_url": "https://www.cancer.gov/publications/dictionaries/cancer-terms/def/mekinist",
          "scope": "Selected documented pharmacological targets; not an exhaustive target profile or PLATO output."
        },
        {
          "source": "drug:Trametinib",
          "target": "gene:MAP2K2",
          "type": "known_target",
          "evidence_status": "literature_annotation",
          "source_url": "https://www.cancer.gov/publications/dictionaries/cancer-terms/def/mekinist",
          "scope": "Selected documented pharmacological targets; not an exhaustive target profile or PLATO output."
        },
        {
          "source": "gene:BRAF",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "P15056",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:MAP2K1",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "Q02750",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:MAP2K2",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "P36507",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:MAPK1",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "P28482",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:MAPK3",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "P27361",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:RAF1",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "P04049",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:KSR1",
          "target": "pathway:R-HSA-5673001",
          "type": "belongs_to_pathway",
          "source_url": "https://reactome.org/content/detail/R-HSA-5673001",
          "reference_accession": "Q8IVT5",
          "membership_basis": "Reactome pathway reference-entity participant",
          "source_query": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities"
        },
        {
          "source": "gene:BRAF",
          "target": "disease_module:literature-context",
          "type": "belongs_to_disease_module",
          "source_url": "https://pubmed.ncbi.nlm.nih.gov/23020132/",
          "membership_basis": "Manually selected BRAF V600 trial inclusion context; single-gene literature background, not a disease-resource gene set."
        }
      ],
      "semantics": "Public reference annotations; undirected interactions, drug-target associations and membership are distinct. Layout positions encode no biology."
    }
  },
  "local_evidence": {
    "cross_drug_paths": [
      {
        "drug": "Dabrafenib",
        "destination_drug": "Trametinib",
        "path": [
          "BRAF",
          "MAP2K1"
        ],
        "length": 1,
        "selection": "editorial existing-edge route"
      },
      {
        "drug": "Dabrafenib",
        "destination_drug": "Trametinib",
        "path": [
          "BRAF",
          "MAP2K2"
        ],
        "length": 1,
        "selection": "editorial existing-edge route"
      },
      {
        "drug": "Dabrafenib",
        "destination_drug": "Trametinib",
        "path": [
          "BRAF",
          "RAF1",
          "MAP2K2"
        ],
        "length": 2,
        "selection": "editorial existing-edge route, not a causal bypass claim"
      }
    ],
    "top_paths": [
      {
        "drug": "Dabrafenib",
        "path": [
          "BRAF"
        ],
        "length": 0,
        "selection": "single-gene literature background membership"
      },
      {
        "drug": "Trametinib",
        "path": [
          "MAP2K1",
          "BRAF"
        ],
        "length": 1,
        "selection": "existing-edge route to literature background member"
      },
      {
        "drug": "Trametinib",
        "path": [
          "MAP2K2",
          "BRAF"
        ],
        "length": 1,
        "selection": "existing-edge route to literature background member"
      }
    ]
  },
  "network_connectivity": {
    "status": "reference_annotation_only_not_computed",
    "note": "Paths are selected from the included annotation graph. No interpreter, disease proximity or cell-line calculation was run."
  },
  "network_statistics": {
    "valid": false,
    "quality_status": "not_evaluated",
    "reason": "reference_example_no_model_or_statistics"
  },
  "reference_provenance": {
    "retrieved_on": "2026-10-21",
    "uniprot_entries": [
      {
        "gene": "BRAF",
        "accession": "P15056",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1990-11-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "2004-07-19",
          "entryVersion": 269,
          "sequenceVersion": 4
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/P15056.json",
        "source_sha256": "f4d1e9e4534be377e91691c09016add1972ecb6f149195320fad2b40260cdc46",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "P15056",
              "intActId": "EBI-365980"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q02750",
              "geneName": "MAP2K1",
              "intActId": "EBI-492564"
            },
            "numberOfExperiments": 65,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P15056",
              "intActId": "EBI-365980"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P36507",
              "geneName": "MAP2K2",
              "intActId": "EBI-1056930"
            },
            "numberOfExperiments": 13,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P15056",
              "intActId": "EBI-365980"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P04049",
              "geneName": "RAF1",
              "intActId": "EBI-365996"
            },
            "numberOfExperiments": 73,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "MAP2K1",
        "accession": "Q02750",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1993-07-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "2007-01-23",
          "entryVersion": 249,
          "sequenceVersion": 2
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/Q02750.json",
        "source_sha256": "cf5f54050e2a6f837987900ff7145b48bf4ca7ab1a4e49d2710889bd9e59ba24",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "Q02750",
              "intActId": "EBI-492564"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P15056",
              "geneName": "BRAF",
              "intActId": "EBI-365980"
            },
            "numberOfExperiments": 65,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "Q02750",
              "intActId": "EBI-492564"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P28482",
              "geneName": "MAPK1",
              "intActId": "EBI-959949"
            },
            "numberOfExperiments": 2,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "Q02750",
              "intActId": "EBI-492564"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P27361",
              "geneName": "MAPK3",
              "intActId": "EBI-73995"
            },
            "numberOfExperiments": 2,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "Q02750",
              "intActId": "EBI-492564"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P04049",
              "geneName": "RAF1",
              "intActId": "EBI-365996"
            },
            "numberOfExperiments": 39,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "MAP2K2",
        "accession": "P36507",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1994-06-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "1994-06-01",
          "entryVersion": 254,
          "sequenceVersion": 1
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/P36507.json",
        "source_sha256": "d316777cec1ef826a0ed062f7fc521d9729c782df917a46066331e00949186de",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "P36507",
              "intActId": "EBI-1056930"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P15056",
              "geneName": "BRAF",
              "intActId": "EBI-365980"
            },
            "numberOfExperiments": 13,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P36507",
              "intActId": "EBI-1056930"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q8IVT5",
              "geneName": "KSR1",
              "intActId": "EBI-486984"
            },
            "numberOfExperiments": 10,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P36507",
              "intActId": "EBI-1056930"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P04049",
              "geneName": "RAF1",
              "intActId": "EBI-365996"
            },
            "numberOfExperiments": 7,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "MAPK1",
        "accession": "P28482",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1992-12-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "2007-01-23",
          "entryVersion": 261,
          "sequenceVersion": 3
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/P28482.json",
        "source_sha256": "5166e9247fd19dd7b92dde5db62e68d67d0db808ee9582c078d33ca64ac3a54b",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "P28482",
              "intActId": "EBI-959949"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q02750",
              "geneName": "MAP2K1",
              "intActId": "EBI-492564"
            },
            "numberOfExperiments": 2,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "MAPK3",
        "accession": "P27361",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1992-08-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "2007-01-23",
          "entryVersion": 258,
          "sequenceVersion": 4
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/P27361.json",
        "source_sha256": "34ab28a7c926756d4461e5486239b127c95382ef52515b927296658386caea1a",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "P27361",
              "intActId": "EBI-73995"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q02750",
              "geneName": "MAP2K1",
              "intActId": "EBI-492564"
            },
            "numberOfExperiments": 2,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "RAF1",
        "accession": "P04049",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "1986-11-01",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "1986-11-01",
          "entryVersion": 271,
          "sequenceVersion": 1
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/P04049.json",
        "source_sha256": "6e8cc60c356b3b0578a3c04946af28ecd31aee42960929870baab327dc6697ca",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "P04049",
              "intActId": "EBI-365996"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P15056",
              "geneName": "BRAF",
              "intActId": "EBI-365980"
            },
            "numberOfExperiments": 73,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P04049",
              "intActId": "EBI-365996"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q8IVT5",
              "geneName": "KSR1",
              "intActId": "EBI-486984"
            },
            "numberOfExperiments": 5,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P04049",
              "intActId": "EBI-365996"
            },
            "interactantTwo": {
              "uniProtKBAccession": "Q02750",
              "geneName": "MAP2K1",
              "intActId": "EBI-492564"
            },
            "numberOfExperiments": 39,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "P04049",
              "intActId": "EBI-365996"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P36507",
              "geneName": "MAP2K2",
              "intActId": "EBI-1056930"
            },
            "numberOfExperiments": 7,
            "organismDiffer": false
          }
        ]
      },
      {
        "gene": "KSR1",
        "accession": "Q8IVT5",
        "taxon_id": 9606,
        "entry_audit": {
          "firstPublicDate": "2004-12-21",
          "lastAnnotationUpdateDate": "2026-09-02",
          "lastSequenceUpdateDate": "2014-09-03",
          "entryVersion": 197,
          "sequenceVersion": 3
        },
        "source_url": "https://rest.uniprot.org/uniprotkb/Q8IVT5.json",
        "source_sha256": "d824c8c065adaa9f2521522d137796d779b5e9cc16ce153948f2142ae9b8b22e",
        "selected_interactions": [
          {
            "interactantOne": {
              "uniProtKBAccession": "Q8IVT5",
              "intActId": "EBI-486984"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P36507",
              "geneName": "MAP2K2",
              "intActId": "EBI-1056930"
            },
            "numberOfExperiments": 10,
            "organismDiffer": false
          },
          {
            "interactantOne": {
              "uniProtKBAccession": "Q8IVT5",
              "intActId": "EBI-486984"
            },
            "interactantTwo": {
              "uniProtKBAccession": "P04049",
              "geneName": "RAF1",
              "intActId": "EBI-365996"
            },
            "numberOfExperiments": 5,
            "organismDiffer": false
          }
        ]
      }
    ],
    "reactome": {
      "pathway": "R-HSA-5673001",
      "source_url": "https://reactome.org/ContentService/data/participants/R-HSA-5673001/referenceEntities",
      "source_sha256": "b550b8c29219156f30a0661420838f0ca6c382ae6498630e8f76b2938a692772",
      "selected_participants": {
        "BRAF": "P15056",
        "MAP2K1": "Q02750",
        "MAP2K2": "P36507",
        "MAPK1": "P28482",
        "MAPK3": "P27361",
        "RAF1": "P04049",
        "KSR1": "Q8IVT5"
      }
    }
  }
};
});
