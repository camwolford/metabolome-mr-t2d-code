"""Count the Open Targets Platform GWAS credible sets that contain each instrumenting SNP, and the traits they map to.

Usage: python3 02_annotate_instruments_open_targets.py <out.tsv> <rsID> [<rsID> ...]

The manuscript used Open Targets Platform release 26.09; the API serves the current release, which is printed.
"""

import csv
import json
import sys
import time
import urllib.request

API = "https://api.platform.opentargets.org/api/v4/graphql"
DISEASE_ONTOLOGIES = ("MONDO", "HP", "Orphanet", "DOID")  # disease or phenotype terms; other terms are measurements
CREDIBLE_SETS = """query($v:String!,$i:Int!){variant(variantId:$v){credibleSets(studyTypes:[gwas], page:{index:$i,size:100}){
  count rows{study{diseases{id}}}}}}"""


def query(text, variables=None):
    body = json.dumps({"query": text, "variables": variables or {}}).encode()
    request = urllib.request.Request(API, body, {"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            response = json.load(urllib.request.urlopen(request, timeout=120))
            if "errors" in response:
                raise RuntimeError(response["errors"])
            return response["data"]
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3)


def variant_ids(rsid):
    """Open Targets variant IDs (chromosome_position_ref_alt, GRCh38) carrying this rsID."""
    hits = query(
        'query($s:String!){search(queryString:$s, entityNames:["variant"], page:{index:0,size:10}){hits{id}}}',
        {"s": rsid},
    )["search"]["hits"]
    ids = []
    for hit in hits:
        variant = query("query($v:String!){variant(variantId:$v){id rsIds}}", {"v": hit["id"]})["variant"]
        if variant and rsid in (variant["rsIds"] or []):
            ids.append(variant["id"])
    return ids


def gwas_credible_sets(variant_id):
    rows, page = [], 0
    while True:
        result = query(CREDIBLE_SETS, {"v": variant_id, "i": page})["variant"]["credibleSets"]
        rows += result["rows"]
        if len(rows) >= result["count"] or not result["rows"]:
            return result["count"], rows
        page += 1


def main(out_path, rsids):
    version = query("{meta{dataVersion{year month}}}")["meta"]["dataVersion"]
    print(f"Open Targets Platform release {version['year']}.{version['month']}")
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rsID", "Open_Targets_variant", "GWAS_credible_sets", "Distinct_traits", "Disease_phenotype_terms"])
        for rsid in rsids:
            ids = variant_ids(rsid)
            if not ids:  # not in the Open Targets variant index
                writer.writerow([rsid, "NA", 0, 0, 0])
            for variant_id in ids:
                count, rows = gwas_credible_sets(variant_id)
                traits = {term["id"] for row in rows for term in (row["study"]["diseases"] or [])}
                diseases = {t for t in traits if t.split("_")[0] in DISEASE_ONTOLOGIES}
                writer.writerow([rsid, variant_id, count, len(traits), len(diseases)])


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
