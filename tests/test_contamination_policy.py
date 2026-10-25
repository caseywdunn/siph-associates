"""Consequential matching and evidence-preservation checks for both tiers."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from contamination_policy import ContaminationPolicy

def test_policy():
    policy = ContaminationPolicy(catalog=[{"catalog_id":"BAC00289","genome_id":"GCF_008364245.1","cluster_members":"GCF_008364245.1,MAGSP0001"}, {"catalog_id":"","genome_id":"MAGSP0001","cluster_members":""}, {"catalog_id":"myco","genome_id":"MAGSP0012"}])
    source = {"sample_id":"Ahuja2026__NA38","taxonomy":"d__Bacteria;g__Pseudomonas_E;s__Pseudomonas_E x"}
    assert policy.annotate(source)["contaminant_taxon"] == "true"
    assert "biological_include" not in source
    assert policy.annotate({"analysis_sample_id":"Ahuja2026__NA31", "sample_id":"NA31", "library_id":"Ahuja2026:NA31"})["biological_include"] == "false"
    assert policy.annotate({"sample_id":"Ahuja2026__NA38","taxonomy":"d__Bacteria;g__Pseudoalteromonas;s__Pseudoalteromonas x"})["biological_include"] == "true"
    assert policy.annotate({"sample_id":"Ahuja2026__NA31","target_id":"myco"})["excluded_library"] == "true"
    assert policy.annotate({"sample_id":"Ahuja2026__NA31","target_id":"myco"})["contaminant_taxon"] == "false"
    assert policy.annotate({"sample_id":"Church2025__TMAG-K5632","reporting_unit":"Nematozoa"})["contaminant_taxon"] == "true"
    assert "moderate quality" in policy.annotate({"sample_id":"Church2025__TMAG-K5633","reporting_unit":"Nematozoa"})["contaminant_basis"]
    assert policy.annotate({"sample_id":"Church2025__TMAG-K5632","reporting_unit":"Copepoda"})["biological_include"] == "true"
    assert policy.annotate({"sample_id":"Church2025__TMAG-K5632","target_id":"BAC00289"})["contaminant_taxon"] == "true"
    assert policy.annotate({"sample_id":"Manko2026__T676-SS3"})["biological_include"] == "false"

if __name__ == "__main__":
    test_policy()
    print("Contamination-policy checks passed")
