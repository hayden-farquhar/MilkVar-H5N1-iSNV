"""
Build the SRA manifest for MilkVar (Project 88).

Queries NCBI Entrez for H5N1 cattle/feline/milk SRA runs from key BioProjects,
extracts metadata, and produces a frozen manifest TSV.

Usage:
    python3 scripts/build_manifest.py
"""

import time
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "hayden.farquhar@icloud.com"
Entrez.tool = "MilkVar_manifest_builder"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_DIR = PROJECT_ROOT / "data" / "manifests"
MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

KEY_BIOPROJECTS = [
    "PRJNA1102327",   # USDA H5Nx immediate releases
    "PRJNA1416438",   # California dairy farm surveillance
    "PRJNA1219588",   # USDA H5N1 D1.1 genotype
    "PRJNA1114404",   # HPAI dairy cattle spillover, Cornell (includes feline)
    "PRJNA1312032",   # Pasteurized milk passive surveillance
    "PRJNA1165336",   # St. Jude retail milk
    "PRJNA1165321",   # St. Jude retail milk (companion)
]

HOST_KEYWORDS = [
    "bos taurus", "cattle", "cow", "bovine", "dairy",
    "felis catus", "cat", "feline",
    "milk", "retail milk", "pasteurized", "pasteurised",
]

EXCLUDE_HOSTS = [
    "chicken", "gallus", "duck", "anas", "goose", "anser",
    "turkey", "meleagris", "wild bird", "gull", "pelican",
    "hawk", "eagle", "owl", "crow", "corvus", "sparrow",
]


def search_sra_by_bioproject(bioproject: str) -> list[str]:
    """Search SRA for all runs in a BioProject."""
    query = f"{bioproject}[BioProject] AND influenza[Organism]"
    handle = Entrez.esearch(db="sra", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"  {bioproject}: {len(ids)} SRA IDs found")
    return ids


def fetch_sra_metadata_batch(id_list: list[str], batch_size: int = 200) -> list[dict]:
    """Fetch SRA metadata in batches and parse XML."""
    records = []
    for i in tqdm(range(0, len(id_list), batch_size), desc="  Fetching"):
        batch = id_list[i:i + batch_size]
        handle = Entrez.efetch(db="sra", id=",".join(batch), rettype="full", retmode="xml")
        xml_text = handle.read()
        handle.close()

        if isinstance(xml_text, bytes):
            xml_text = xml_text.decode("utf-8")

        root = ET.fromstring(f"<root>{xml_text}</root>")
        for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
            rec = parse_experiment_package(exp_pkg)
            if rec:
                records.append(rec)

        time.sleep(0.4)  # rate limit

    return records


def parse_experiment_package(pkg) -> dict | None:
    """Parse a single SRA EXPERIMENT_PACKAGE XML element."""
    run_set = pkg.find(".//RUN_SET/RUN")
    if run_set is None:
        return None

    run_accession = run_set.get("accession", "")
    total_bases = run_set.get("total_bases", "")
    total_spots = run_set.get("total_spots", "")

    experiment = pkg.find(".//EXPERIMENT")
    exp_accession = experiment.get("accession", "") if experiment is not None else ""

    sample = pkg.find(".//SAMPLE")
    biosample = ""
    if sample is not None:
        for xref in sample.findall(".//EXTERNAL_ID"):
            if xref.get("namespace") == "BioSample":
                biosample = xref.text or ""
                break

    # Extract host, collection date, geographic location from sample attributes
    host = ""
    collection_date = ""
    geo_loc = ""
    isolate = ""
    sample_attrs = pkg.findall(".//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
    for attr in sample_attrs:
        tag = (attr.find("TAG").text or "").lower() if attr.find("TAG") is not None else ""
        val = (attr.find("VALUE").text or "") if attr.find("VALUE") is not None else ""
        if tag in ("host", "specific_host", "host scientific name"):
            host = val
        elif tag in ("collection_date", "collection date"):
            collection_date = val
        elif tag in ("geo_loc_name", "geographic location"):
            geo_loc = val
        elif tag in ("isolate", "strain"):
            isolate = val

    # Library strategy
    lib_descriptor = pkg.find(".//LIBRARY_DESCRIPTOR")
    library_strategy = ""
    library_source = ""
    if lib_descriptor is not None:
        ls = lib_descriptor.find("LIBRARY_STRATEGY")
        library_strategy = ls.text if ls is not None else ""
        lsrc = lib_descriptor.find("LIBRARY_SOURCE")
        library_source = lsrc.text if lsrc is not None else ""

    # Platform
    platform_elem = pkg.find(".//PLATFORM")
    platform = ""
    instrument = ""
    if platform_elem is not None:
        for child in platform_elem:
            platform = child.tag
            inst = child.find("INSTRUMENT_MODEL")
            instrument = inst.text if inst is not None else ""
            break

    # BioProject
    study = pkg.find(".//STUDY")
    bioproject = ""
    if study is not None:
        for xref in study.findall(".//EXTERNAL_ID"):
            if xref.get("namespace") == "BioProject":
                bioproject = xref.text or ""
                break

    return {
        "run_accession": run_accession,
        "experiment_accession": exp_accession,
        "biosample": biosample,
        "bioproject": bioproject,
        "host": host,
        "collection_date": collection_date,
        "geo_loc_name": geo_loc,
        "isolate": isolate,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "platform": platform,
        "instrument": instrument,
        "total_bases": total_bases,
        "total_spots": total_spots,
    }


def classify_host(host_str: str) -> str:
    """Classify host string into categories."""
    h = host_str.lower()
    if any(k in h for k in ["bos taurus", "cattle", "cow", "bovine", "dairy"]):
        return "cattle"
    if any(k in h for k in ["felis catus", "cat", "feline"]):
        return "feline"
    if any(k in h for k in ["milk", "retail", "pasteurized", "pasteurised"]):
        return "retail_milk"
    return "other"


def is_target_host(host_str: str) -> bool:
    """Check if the host is a target (cattle/feline/milk) and not excluded."""
    h = host_str.lower()
    if any(ex in h for ex in EXCLUDE_HOSTS):
        return False
    return any(k in h for k in HOST_KEYWORDS)


def main():
    print("=" * 60)
    print("MilkVar Manifest Builder")
    print(f"Date: {date.today().isoformat()}")
    print("=" * 60)

    all_records = []

    for bp in KEY_BIOPROJECTS:
        print(f"\nQuerying {bp}...")
        ids = search_sra_by_bioproject(bp)
        if ids:
            records = fetch_sra_metadata_batch(ids)
            all_records.extend(records)
            print(f"  Parsed {len(records)} experiment packages")

    print(f"\nTotal raw records: {len(all_records)}")

    df = pd.DataFrame(all_records)
    if df.empty:
        print("ERROR: No records retrieved. Check network/Entrez access.")
        return

    # Deduplicate by run_accession
    df = df.drop_duplicates(subset="run_accession")
    print(f"After dedup: {len(df)} unique runs")

    # Classify hosts
    df["host_category"] = df["host"].apply(classify_host)
    df["is_target"] = df["host"].apply(is_target_host)

    # Filter to target hosts
    target_df = df[df["is_target"]].copy()
    print(f"Target-host runs (cattle/feline/milk): {len(target_df)}")

    # Summary
    print("\nHost category breakdown:")
    print(target_df["host_category"].value_counts().to_string())
    print("\nBioProject breakdown:")
    print(target_df["bioproject"].value_counts().to_string())
    print("\nPlatform breakdown:")
    print(target_df["platform"].value_counts().to_string())
    print("\nLibrary strategy breakdown:")
    print(target_df["library_strategy"].value_counts().to_string())

    # Estimate storage
    target_df["total_bases"] = pd.to_numeric(target_df["total_bases"], errors="coerce")
    total_gb = target_df["total_bases"].sum() / 1e9
    print(f"\nEstimated total bases: {total_gb:.1f} Gb")
    print(f"Estimated compressed FASTQ: ~{total_gb * 0.3:.1f} GB (assuming ~3:1 compression)")

    # Save manifest
    today = date.today().strftime("%Y%m%d")
    manifest_path = MANIFEST_DIR / f"v1_{today}.tsv"
    target_df.to_csv(manifest_path, sep="\t", index=False)
    print(f"\nManifest saved: {manifest_path}")
    print(f"Manifest size: {len(target_df)} runs")

    # Also save the full (unfiltered) manifest for reference
    full_path = MANIFEST_DIR / f"v1_{today}_all_hosts.tsv"
    df.to_csv(full_path, sep="\t", index=False)
    print(f"Full manifest (all hosts): {full_path} ({len(df)} runs)")


if __name__ == "__main__":
    main()
