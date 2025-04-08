# python

import csv
import gzip
import json
import logging
import os
import re
import shutil
from json import JSONDecodeError
from typing import Any

import numpy.typing as npt
import pandas as pd
from tqdm import tqdm

from aizynth.classifier.classifier import get_fp, calc_jaccard_similarity

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Add console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def load_smirks(path: str) -> pd.DataFrame:
    logger.info(f"Loading smirks from {path}")
    try:
        with open(path, "r") as f:
            smirks = json.load(f)
            smirks_list = [{"name": name, "smirks": reaction}
                           for name, reaction in smirks.items()]
            smirks = pd.DataFrame(smirks_list)
    except JSONDecodeError:
        smirks = pd.read_json(path, lines=True)
    return smirks


def normalize_name(name: str) -> str:
    """Normalize reaction name for better matching"""
    # Convert to lowercase
    name = name.lower()
    # Remove special characters and replace with spaces
    name = re.sub(r'[^a-z0-9]', ' ', name)
    # Replace multiple spaces with a single space
    name = re.sub(r'\s+', ' ', name)
    # Remove common words that don't help with matching
    stopwords = ['reaction', 'synthesis', 'with', 'and', 'or', 'via', 'from', 'to', 'the', 'of', 'by']
    for word in stopwords:
        name = re.sub(r'\b' + word + r'\b', '', name)
    # Strip leading/trailing spaces and normalize internal spaces
    name = re.sub(r'\s+', ' ', name).strip()
    return name


def calculate_name_similarity(name1: str, name2: str) -> float:
    """Calculate similarity between two reaction names with emphasis on substring matching"""
    # Normalize names
    norm1 = normalize_name(name1)
    norm2 = normalize_name(name2)

    # Check for direct substring matching first
    if norm1 in norm2:
        # Proportional to how much of norm2 is covered by norm1
        return len(norm1) / len(norm2)
    elif norm2 in norm1:
        # Proportional to how much of norm1 is covered by norm2
        return len(norm2) / len(norm1)

    # Fall back to simple token overlap for cases where neither is a substring
    tokens1 = set(token for token in norm1.split() if len(token) > 2)  # filter out very short
    tokens2 = set(token for token in norm2.split() if len(token) > 2)

    if not tokens1 or not tokens2:
        return 0.0

    # Simple token overlap
    intersection = tokens1.intersection(tokens2)
    return len(intersection) / max(len(tokens1), len(tokens2))


def calculate_rxn_similarity(reaction1: str, reaction2: str) -> float:
    """Calculate similarity between two reactions."""
    fp1: npt.NDArray[Any] = get_fp(reaction1, fp="MACCS", concatenate=True)
    fp2: npt.NDArray[Any] = get_fp(reaction2, fp="MACCS", concatenate=True)
    sim: float = calc_jaccard_similarity(fp1, fp2)
    return sim


def group_smirks(labels: pd.DataFrame, to_group: pd.DataFrame, how: str = 'by_score') -> tuple[list[dict], list[dict]]:
    """Group the to_group smirks using the labels smirks. Return grouped and not_grouped smirks"""
    grouped: list = []
    matched_indices = set()

    if how == "by_score":
        not_grouped_indices = set()
        highest_scoring: dict = {}
        # score each to_group reaction similarity to all label smirks, then group after using highest avg sim
        scored: list[dict] = []
        for i, tgrow in tqdm(to_group.iterrows(), total=len(to_group)):
            tgname = tgrow['name']
            tgsmirks = tgrow['smirks']
            scores: dict = {}
            scores['name'] = tgname
            scores['smirks'] = tgsmirks
            label_scored: list[dict] = []
            max_score_rxn: dict = {}
            for _, row in labels.iterrows():
                this_label: dict = {}
                label_name = row['name']
                label_smirks = row['smirks']
                this_label['name'] = label_name
                this_label['smirks'] = label_smirks
                this_label['rxn_sim'] = calculate_rxn_similarity(label_smirks, tgsmirks)
                this_label['name_sim'] = calculate_name_similarity(label_name, tgname)
                this_label['avg_sim'] = (this_label['rxn_sim'] + this_label['name_sim']) / 2
                label_scored.append(this_label)
                if this_label['avg_sim'] >= 0.4:  # baseline threshold
                    if not max_score_rxn or this_label['rxn_sim'] > max_score_rxn['rxn_sim']:
                        max_score_rxn = this_label
            if max_score_rxn:
                # add to highest scoring
                tginfo = {
                    'name': tgname,
                    'smirks': tgsmirks,
                    'rxn_sim': max_score_rxn['rxn_sim'],
                    'name_sim': max_score_rxn['name_sim'],
                    'avg_sim': max_score_rxn['avg_sim'],
                }
            # Check if key exists first
            if max_score_rxn and max_score_rxn['name'] not in highest_scoring:
                highest_scoring[max_score_rxn['name']] = [tginfo]
            elif max_score_rxn and max_score_rxn['name'] in highest_scoring:
                highest_scoring[max_score_rxn['name']].append(tginfo)
            else:
                not_grouped_indices.add(i)
            scores['compared'] = label_scored
            scored.append(scores)
        save_to_json(scored, '../data/dev/scored_smirks.json')

        for _, row in labels.iterrows():
            labelled_rxn: dict = {}
            labelled_rxn['name'] = row['name']
            labelled_rxn['smirks'] = row['smirks']
            if row['name'] in highest_scoring:
                labelled_rxn['matched']: list[dict] = highest_scoring[row['name']]
            else:
                labelled_rxn['matched']: list[dict] = []
            grouped.append(labelled_rxn)
        not_grouped = [dict(to_group.iloc[i]) for i in not_grouped_indices]
        all_matched = sum(len(g['matched']) for g in grouped)
        logger.info(f"Grouped {all_matched} reactions (with repeats)")
        logger.info(f"Not grouped {len(not_grouped)} reactions")

        return grouped, not_grouped

    else:  # by similarity thresholding
        for _, row in labels.iterrows():
            labelled_rxn: dict = {}
            labelled_rxn['name'] = row['name']
            labelled_rxn['smirks'] = row['smirks']
            label_name = row['name']
            label_smirks = row['smirks']
            found: list[dict] = []
            for i, tgrow in to_group.iterrows():
                tgname = tgrow['name']
                tgsmirks = tgrow['smirks']
                similarity: float = calculate_name_similarity(label_name, tgname)
                rxn_similarity: float = calculate_rxn_similarity(label_smirks, tgsmirks)
                if similarity >= 0.4 or rxn_similarity >= 0.5:
                    found.append({"name": tgname, "smirks": tgsmirks})
                    matched_indices.add(i)
            labelled_rxn['matched'] = found
            grouped.append(labelled_rxn)
        not_grouped = [dict(to_group.iloc[i]) for i in range(len(to_group)) if i not in matched_indices]
        # Log summary
        all_matched = sum(len(g['matched']) for g in grouped)
        logger.info(f"Grouped {all_matched} reactions (with repeats)")
        logger.info(f"Not grouped {len(not_grouped)} reactions")
    return grouped, not_grouped


def save_to_json(data: list[dict], path: str) -> bool:
    try:
        with open(path, "w") as f:
            json.dump(data, f, indent=4)
            return True
    except Exception as e:
        logger.error(f"Failed to save data to json: {e}")
        return False


def load_template_codes(home_dir: str) -> dict:
    """Load templates."""
    if not os.path.exists(os.path.join(home_dir, 'uspto_templates.csv')) and os.path.exists(
            os.path.join(home_dir, 'uspto_templates.csv.gz')):
        with gzip.open(os.path.join(home_dir, 'uspto_templates.csv.gz'), 'rb') as f_in:
            with open(os.path.join(home_dir, 'uspto_templates.csv'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Unzipped {home_dir}/uspto_templates.csv.gz to {home_dir}/uspto_templates.csv")
    elif not os.path.exists(os.path.join(home_dir, 'uspto_templates.csv')) and not os.path.exists(
            os.path.join(home_dir, 'uspto_templates.csv.gz')):
        logger.error(
            f"Missing template file! Should be at {home_dir}/uspto_templates.csv or {home_dir}/uspto_templates.csv.gz")

    with open(os.path.join(home_dir, 'uspto_templates.csv'), 'r', encoding='utf-8') as f:
        try:
            template_dict = {}
            reader = csv.DictReader(f, delimiter='\t')
            # Create dictionary from the CSV
            for row in reader:
                if 'template_code' in row and 'retro_template' in row:
                    template_dict[int(row['template_code'])] = row['retro_template']
            return template_dict
        except UnicodeDecodeError:
            logger.error("Failed to decode uspto_templates.csv. Please check the file.")
            return {}


def transform_retro_to_forward(reaction: str) -> str:
    """Given a retrosynthesis reaction, transform it into a forward one."""
    parts = reaction.split('>>')
    product, reactants = parts
    return f"{reactants}>>{product}"


def find_matches(smirks: str, codes: dict) -> list[dict]:
    """Find matches in the template codes."""
    matches: list[dict] = []
    for code, template in tqdm(codes.items()):
        rxn_sim = calculate_rxn_similarity(smirks, transform_retro_to_forward(template))
        if rxn_sim >= 0.46:  # catch as many as possible
            matches.append({"template_code": code, "template": template, "rxn_sim": rxn_sim})
    # Sort matches by similarity
    matches.sort(key=lambda x: x['rxn_sim'], reverse=True)
    logger.info(f"Found {len(matches)} matches for smirks: {smirks}")
    return matches


def match_template_codes(grouped: list[dict], codes: dict) -> None:
    """Match template codes to grouped smirks."""
    logger.info("Matching template codes to grouped smirks")
    for group in tqdm(grouped):
        if group['name'] != 'Amidation':
            continue
        logger.info(f"Matching {group['name']}")
        templates: list[dict] = find_matches(group['smirks'], codes)
        group['uspto_templates'] = templates
        for rxn in group['matched']:
            smirks = rxn['smirks']
            templates: list[dict] = find_matches(smirks, codes)
            rxn['uspto_templates'] = templates
    logger.info("Finished matching template codes")
    save_to_json(grouped, '../data/dev/grouped_smirks_with_codes.json')


def main():
    if not os.path.exists("../data/dev/grouped_smirks.json"):
        syn: pd.DataFrame = load_smirks('../../syndirella/constants/RXN_SMIRKS_CONSTANTS.json')
        rxni: pd.DataFrame = load_smirks('../data/dev/smirks.json')
        grouped, not_grouped = group_smirks(syn, rxni)
        if not save_to_json(grouped, '../data/dev/grouped_smirks.json'):
            logger.error("Failed to save grouped smirks")
        if not save_to_json(not_grouped, '../data/dev/not_grouped_smirks.json'):
            logger.error("Failed to save not_grouped smirks")
    else:
        with open(os.path.join("../data/dev/grouped_smirks.json"), "r") as f:
            grouped = json.load(f)
        codes: dict = load_template_codes(home_dir='..')
        if codes:
            match_template_codes(grouped, codes)


if __name__ == '__main__':
    main()
